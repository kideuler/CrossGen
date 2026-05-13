// ViewerMain.cxx – Qt6 application entry point for the CrossGen viewer.

#include <QApplication>
#include <QSurfaceFormat>
#include <QMessageBox>

#include <iostream>
#include <stdexcept>

#include "viewer/CrossGenWidget.hxx"

int main(int argc, char **argv) {
    // Request a legacy (compatibility) OpenGL profile so that all the
    // existing fixed-function GL draw calls (glBegin/glEnd, glOrtho, …) keep working.
    QSurfaceFormat fmt;
    fmt.setProfile(QSurfaceFormat::CompatibilityProfile);
    fmt.setVersion(2, 1);
    fmt.setSamples(8); // multisampling (MSAA)
    QSurfaceFormat::setDefaultFormat(fmt);

    QApplication app(argc, argv);
    app.setApplicationName("CrossGen Viewer");

    if (argc < 2) {
        std::cerr << "Usage: Viewer <mesh.obj>\n";
        QMessageBox::critical(nullptr, "CrossGen Viewer",
                              "Usage: Viewer <mesh.obj>\nNo mesh file provided.");
        return 1;
    }

    std::string meshPath = argv[1];

    CrossGenWidget *widget = nullptr;
    try {
        widget = new CrossGenWidget(meshPath);
    } catch (const std::exception &e) {
        std::cerr << "Failed to load mesh: " << e.what() << "\n";
        QMessageBox::critical(nullptr, "CrossGen Viewer",
                              QString("Failed to load mesh:\n%1").arg(e.what()));
        return 2;
    }

    widget->setWindowTitle(QString("CrossGen Viewer – %1")
                           .arg(QString::fromStdString(meshPath)));
    widget->resize(1200, 900);
    widget->show();

    return app.exec();
}
