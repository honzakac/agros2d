// This file is part of Agros2D.
//
// Agros2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Agros2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Agros2D.  If not, see <http://www.gnu.org/licenses/>.
//
// hp-FEM group (http://hpfem.org/)
// University of Nevada, Reno (UNR) and University of West Bohemia, Pilsen
// Email: agros2d@googlegroups.com, home page: http://hpfem.org/agros2d/

#include "config.h"

Config::Config()
{
    load();
}

Config::~Config()
{
    save();
}

void Config::load()
{
    QSettings settings;

    // general
    guiStyle = settings.value("General/GUIStyle").toString();
    language = settings.value("General/Language", QLocale::system().name()).toString();
    defaultPhysicField = (PhysicField) settings.value("General/DefaultPhysicField", PhysicField_Electrostatic).toInt();
    if (defaultPhysicField == PhysicField_Undefined) defaultPhysicField = PhysicField_Electrostatic;

    checkVersion = settings.value("General/CheckVersion", true).toBool();
    showConvergenceChart = settings.value("General/ShowConvergenceChart", true).toBool();
    enabledApplicationLog = settings.value("General/EnabledApplicationLog", false).toBool();
    enabledProgressLog = settings.value("General/EnabledProgressLog", true).toBool();
    lineEditValueShowResult = settings.value("General/LineEditValueShowResult", false).toBool();
    saveProblemWithSolution = settings.value("Solver/SaveProblemWithSolution", false).toBool();

    // zoom
    zoomToMouse = settings.value("Geometry/ZoomToMouse", true).toBool();

    // delete files
    deleteTriangleMeshFiles = settings.value("Solver/DeleteTriangleMeshFiles", true).toBool();
    deleteHermes2DMeshFile = settings.value("Solver/DeleteHermes2DMeshFile", true).toBool();

    // colors
    colorBackground = settings.value("SceneViewSettings/ColorBackground", COLORBACKGROUND).value<QColor>();
    colorGrid = settings.value("SceneViewSettings/ColorGrid", COLORGRID).value<QColor>();
    colorCross = settings.value("SceneViewSettings/ColorCross", COLORCROSS).value<QColor>();
    colorNodes = settings.value("SceneViewSettings/ColorNodes", COLORNODES).value<QColor>();
    colorEdges = settings.value("SceneViewSettings/ColorEdges", COLOREDGES).value<QColor>();
    colorLabels = settings.value("SceneViewSettings/ColorLabels", COLORLABELS).value<QColor>();
    colorContours = settings.value("SceneViewSettings/ColorContours", COLORCONTOURS).value<QColor>();
    colorVectors = settings.value("SceneViewSettings/ColorVectors", COLORVECTORS).value<QColor>();
    colorInitialMesh = settings.value("SceneViewSettings/ColorInitialMesh", COLORINITIALMESH).value<QColor>();
    colorSolutionMesh = settings.value("SceneViewSettings/ColorSolutionMesh", COLORSOLUTIONMESH).value<QColor>();
    colorHighlighted = settings.value("SceneViewSettings/ColorHighlighted", COLORHIGHLIGHTED).value<QColor>();
    colorSelected = settings.value("SceneViewSettings/ColorSelected", COLORSELECTED).value<QColor>();

    // geometry
    nodeSize = settings.value("SceneViewSettings/NodeSize", 6.0).toDouble();
    edgeWidth = settings.value("SceneViewSettings/EdgeWidth", 2.0).toDouble();
    labelSize = settings.value("SceneViewSettings/LabelSize", 6.0).toDouble();
    angleSegmentsCount = settings.value("Geometry/AngleSegmentsCount", 5).toInt();

    // grid
    gridStep = settings.value("SceneViewSettings/GridStep", 0.05).toDouble();
    // rulers
    showRulers = settings.value("SceneViewSettings/ShowRulers", false).toBool();
    // snap to grid
    snapToGrid = settings.value("SceneViewSettings/SnapToGrid", false).toBool();

    // countour
    contoursCount = settings.value("SceneViewSettings/ContoursCount", 15).toInt();

    // scalar view
    paletteType = (PaletteType) settings.value("SceneViewSettings/PaletteType", Palette_Jet).toInt();
    paletteFilter = settings.value("SceneViewSettings/PaletteFilter", false).toBool();
    paletteSteps = settings.value("SceneViewSettings/PaletteSteps", 30).toInt();
    scalarRangeLog = settings.value("SceneViewSettings/ScalarRangeLog", false).toBool();
    scalarRangeBase = settings.value("SceneViewSettings/ScalarRangeBase", 10).toDouble();

    // vector view
    vectorProportional = settings.value("SceneViewSettings/VectorProportional", false).toBool();
    vectorColor = settings.value("SceneViewSettings/VectorColor", true).toBool();
    vectorCount = settings.value("SceneViewSettings/VectorNumber", 50).toInt();
    vectorScale = settings.value("SceneViewSettings/VectorScale", 0.6).toDouble();

    // order view
    orderLabel = settings.value("SceneViewSettings/OrderLabel", false).toBool();

    // 3d
    scalarView3DLighting = settings.value("SceneViewSettings/ScalarView3DLighting", false).toBool();
    scalarView3DAngle = settings.value("SceneViewSettings/ScalarView3DAngle", 270).toDouble();
    scalarView3DBackground = settings.value("SceneViewSettings/ScalarView3DBackground", true).toBool();

    // adaptivity
    isoOnly = settings.value("Adaptivity/IsoOnly", ADAPTIVITY_ISOONLY).toBool();
    convExp = settings.value("Adaptivity/ConvExp", ADAPTIVITY_CONVEXP).toDouble();
    threshold = settings.value("Adaptivity/Threshold", ADAPTIVITY_THRESHOLD).toDouble();
    strategy = settings.value("Adaptivity/Strategy", ADAPTIVITY_STRATEGY).toInt();
    meshRegularity = settings.value("Adaptivity/MeshRegularity", ADAPTIVITY_MESHREGULARITY).toInt();

    // command argument
    commandTriangle = settings.value("Commands/Triangle", COMMANDS_TRIANGLE).toString();
    commandFFmpeg = settings.value("Commands/FFmpeg", COMMANDS_FFMPEG).toString();

    // global script
    globalScript = settings.value("Python/GlobalScript", "").toString();
}

void Config::save()
{
    QSettings settings;

    // general
    settings.setValue("General/GUIStyle", guiStyle);
    settings.setValue("General/Language", language);
    settings.setValue("General/DefaultPhysicField", defaultPhysicField);

    settings.setValue("General/CheckVersion", checkVersion);
    settings.setValue("General/ShowConvergenceChart", showConvergenceChart);
    settings.setValue("General/EnabledApplicationLog", enabledApplicationLog);
    settings.setValue("General/EnabledProgressLog", enabledProgressLog);
    settings.setValue("General/LineEditValueShowResult", lineEditValueShowResult);
    settings.setValue("General/SaveProblemWithSolution", saveProblemWithSolution);

    // zoom
    settings.setValue("General/ZoomToMouse", zoomToMouse);

    // delete files
    settings.setValue("Solver/DeleteTriangleMeshFiles", deleteTriangleMeshFiles);
    settings.setValue("Solver/DeleteHermes2DMeshFile", deleteHermes2DMeshFile);

    // colors
    settings.setValue("SceneViewSettings/ColorBackground", colorBackground);
    settings.setValue("SceneViewSettings/ColorGrid", colorGrid);
    settings.setValue("SceneViewSettings/ColorCross", colorCross);
    settings.setValue("SceneViewSettings/ColorNodes", colorNodes);
    settings.setValue("SceneViewSettings/ColorEdges", colorEdges);
    settings.setValue("SceneViewSettings/ColorLabels", colorLabels);
    settings.setValue("SceneViewSettings/ColorContours", colorContours);
    settings.setValue("SceneViewSettings/ColorVectors", colorVectors);
    settings.setValue("SceneViewSettings/ColorInitialMesh", colorInitialMesh);
    settings.setValue("SceneViewSettings/ColorSolutionMesh", colorSolutionMesh);
    settings.setValue("SceneViewSettings/ColorInitialMesh", colorHighlighted);
    settings.setValue("SceneViewSettings/ColorSolutionMesh", colorSelected);

    // geometry
    settings.setValue("SceneViewSettings/NodeSize", nodeSize);
    settings.setValue("SceneViewSettings/EdgeWidth", edgeWidth);
    settings.setValue("SceneViewSettings/LabelSize", labelSize);
    settings.setValue("SceneViewSettings/AngleSegmentsCount", angleSegmentsCount);

    // grid
    settings.setValue("SceneViewSettings/GridStep", gridStep);
    // rulers
    settings.setValue("SceneViewSettings/ShowRulers", showRulers);
    // snap to grid
    settings.setValue("SceneViewSettings/SnapToGrid", snapToGrid);

    // countour
    settings.setValue("SceneViewSettings/ContoursCount", contoursCount);

    // scalar view
    settings.setValue("SceneViewSettings/PaletteType", paletteType);
    settings.setValue("SceneViewSettings/PaletteFilter", paletteFilter);
    settings.setValue("SceneViewSettings/PaletteSteps", paletteSteps);
    settings.setValue("SceneViewSettings/ScalarRangeLog", scalarRangeLog);
    settings.setValue("SceneViewSettings/ScalarRangeBase", scalarRangeBase);

    // vector view
    settings.setValue("SceneViewSettings/VectorProportional", vectorProportional);
    settings.setValue("SceneViewSettings/VectorColor", vectorColor);
    settings.setValue("SceneViewSettings/VectorNumber", vectorCount);
    settings.setValue("SceneViewSettings/VectorScale", vectorScale);

    // order view
    settings.setValue("SceneViewSettings/OrderLabel", orderLabel);

    // 3d
    settings.setValue("SceneViewSettings/ScalarView3DLighting", scalarView3DLighting);
    settings.setValue("SceneViewSettings/ScalarView3DAngle", scalarView3DAngle);
    settings.setValue("SceneViewSettings/ScalarView3DBackground", scalarView3DBackground);

    // adaptivity
    settings.setValue("Adaptivity/IsoOnly", isoOnly);
    settings.setValue("Adaptivity/ConvExp", convExp);
    settings.setValue("Adaptivity/Threshold", threshold);
    settings.setValue("Adaptivity/Strategy", strategy);
    settings.setValue("Adaptivity/MeshRegularity", meshRegularity);

    // command argument
    settings.setValue("Commands/Triangle", commandTriangle);
    settings.setValue("Commands/FFmpeg", commandFFmpeg);

    // global script
    settings.setValue("Python/GlobalScript", globalScript);
}
