#include "UVIsoTrace.hxx"

UVIsoTrace::UVIsoTrace(std::shared_ptr<UVGParam> uvParam, int nU, int nV)
    : uvParam_(uvParam)
    , uIntervalTree_(*uvParam, IntervalTree::Axis::U)
    , vIntervalTree_(*uvParam, IntervalTree::Axis::V)
{
    nU_ = nU;
    nV_ = nV;
    const Eigen::VectorXd &u = uvParam_->getU();
    const Eigen::VectorXd &v = uvParam_->getV();
    if (u.size() == 0 || v.size() == 0) {
        uvMin_ = {0.0, 0.0};
        uvMax_ = {1.0, 1.0};
    } else {
        uvMin_ = {u.minCoeff(), v.minCoeff()};
        uvMax_ = {u.maxCoeff(), v.maxCoeff()};
    }

    // Compute step sizes for uniform sampling
    deltaU_ = (uvMax_[0] - uvMin_[0]) / static_cast<double>(nU_-1);
    deltaV_ = (uvMax_[1] - uvMin_[1]) / static_cast<double>(nV_-1);

    // Compute which triangles are valid for tracing (not near singularities or flipped)
    const Mesh& mesh = uvParam_->getCutMesh().getCutMesh();
    int nT = static_cast<int>(mesh.triangles.size());
    isValidTriangle_.resize(nT, true);
    for (int t = 0; t < nT; ++t) {
        const Triangle& tri = mesh.triangles[t];
        int i = tri[0], j = tri[1], k = tri[2];

        // Check if triangle is flipped in UV space
        double signedArea2 = (u(j) - u(i)) * (v(k) - v(i)) - (u(k) - u(i)) * (v(j) - v(i));
        if (signedArea2 <= 0.0) {
            isValidTriangle_[t] = false;
            continue;
        }

        // check if any vertex within this triangle is singular
        const std::vector<bool>& isSingularVertex = uvParam_->getCutMesh().getIsSingularVertex();
        if (isSingularVertex[i] || isSingularVertex[j] || isSingularVertex[k]) {
            isValidTriangle_[t] = false;
            continue;
        }
    }
}


void UVIsoTrace::printQueries() const {

    //print u queries
    for (double uVal = uvMin_[0]; uVal <= uvMax_[0]; uVal += deltaU_) {
        std::vector<int> triIndices = uIntervalTree_.query(uVal);
        std::cout << "u = " << uVal << ": " << triIndices.size() << " intersecting triangles\n";
        // print the triangle indices for debugging
        for (int idx : triIndices) {
            std::cout << "  Triangle " << idx << (isValidTriangle_[idx] ? " (valid)" : " (invalid)") << "\n";
        }   
    }

    // print v queries
    for (double vVal = uvMin_[1]; vVal <= uvMax_[1]; vVal += deltaV_) {
        std::vector<int> triIndices = vIntervalTree_.query(vVal);
        std::cout << "v = " << vVal << ": " << triIndices.size() << " intersecting triangles\n";
        // print the triangle indices for debugging
        for (int idx : triIndices) {
            std::cout << "  Triangle " << idx << (isValidTriangle_[idx] ? " (valid)" : " (invalid)") << "\n";
        }      
    }
}