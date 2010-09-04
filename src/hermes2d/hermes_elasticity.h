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

#ifndef HERMES_ELASTICITY_H
#define HERMES_ELASTICITY_H

#include "util.h"
#include "hermes_field.h"

struct HermesElasticity : public HermesField
{
    Q_OBJECT
public:
    HermesElasticity() { physicField = PhysicField_Elasticity; }

    inline int numberOfSolution() { return 2; }
    bool hasHarmonic() { return false; }
    bool hasTransient() { return false; }

    void readEdgeMarkerFromDomElement(QDomElement *element);
    void writeEdgeMarkerToDomElement(QDomElement *element, SceneEdgeMarker *marker);
    void readLabelMarkerFromDomElement(QDomElement *element);
    void writeLabelMarkerToDomElement(QDomElement *element, SceneLabelMarker *marker);

    LocalPointValue *localPointValue(Point point);
    QStringList localPointValueHeader();

    SurfaceIntegralValue *surfaceIntegralValue();
    QStringList surfaceIntegralValueHeader();

    VolumeIntegralValue *volumeIntegralValue();
    QStringList volumeIntegralValueHeader();

    inline bool physicFieldBCCheck(PhysicFieldBC physicFieldBC) { return (physicFieldBC == PhysicFieldBC_Elasticity_Fixed ||
                                                                          physicFieldBC == PhysicFieldBC_Elasticity_Free); }
    inline bool physicFieldVariableCheck(PhysicFieldVariable physicFieldVariable) { return (physicFieldVariable == PhysicFieldVariable_Elasticity_VonMisesStress); }

    SceneEdgeMarker *newEdgeMarker();
    SceneEdgeMarker *newEdgeMarker(PyObject *self, PyObject *args);
    SceneLabelMarker *newLabelMarker();
    SceneLabelMarker *newLabelMarker(PyObject *self, PyObject *args);

    QList<SolutionArray *> *solve(ProgressItemSolve *progressItemSolve);

    inline PhysicFieldVariable contourPhysicFieldVariable() { return PhysicFieldVariable_Elasticity_VonMisesStress; }
    inline PhysicFieldVariable scalarPhysicFieldVariable() { return PhysicFieldVariable_Elasticity_VonMisesStress; }
    inline PhysicFieldVariableComp scalarPhysicFieldVariableComp() { return PhysicFieldVariableComp_Scalar; }
    inline PhysicFieldVariable vectorPhysicFieldVariable() { return PhysicFieldVariable_Undefined; }

    void fillComboBoxScalarVariable(QComboBox *cmbFieldVariable)
    {
        cmbFieldVariable->addItem(physicFieldVariableString(PhysicFieldVariable_Elasticity_VonMisesStress), PhysicFieldVariable_Elasticity_VonMisesStress);
    }

    void fillComboBoxVectorVariable(QComboBox *cmbFieldVariable)
    {
        // FIX ME
        cmbFieldVariable->addItem(physicFieldVariableString(PhysicFieldVariable_Heat_TemperatureGradient), PhysicFieldVariable_Heat_TemperatureGradient);
    }

    void showLocalValue(QTreeWidget *trvWidget, LocalPointValue *localPointValue);
    void showSurfaceIntegralValue(QTreeWidget *trvWidget, SurfaceIntegralValue *surfaceIntegralValue);
    void showVolumeIntegralValue(QTreeWidget *trvWidget, VolumeIntegralValue *volumeIntegralValue);

    ViewScalarFilter *viewScalarFilter(PhysicFieldVariable physicFieldVariable, PhysicFieldVariableComp physicFieldVariableComp);
};

class LocalPointValueElasticity : public LocalPointValue
{
public:
    double young_modulus;
    double poisson_ratio;
    double von_mises_stress;

    LocalPointValueElasticity(Point &point);
    double variableValue(PhysicFieldVariable physicFieldVariable, PhysicFieldVariableComp physicFieldVariableComp);
    QStringList variables();
};

class SurfaceIntegralValueElasticity : public SurfaceIntegralValue
{
protected:
    void calculateVariables(int i);

public:
    SurfaceIntegralValueElasticity();

    QStringList variables();
};

class VolumeIntegralValueElasticity : public VolumeIntegralValue
{
protected:
    void calculateVariables(int i);
    void initSolutions();

public:
    VolumeIntegralValueElasticity();
    QStringList variables();
};

class ViewScalarFilterElasticity : public ViewScalarFilter
{
public:
    ViewScalarFilterElasticity(Tuple<MeshFunction *> sln, PhysicFieldVariable physicFieldVariable, PhysicFieldVariableComp physicFieldVariableComp) :
            ViewScalarFilter(sln, physicFieldVariable, physicFieldVariableComp) {}

protected:
    void calculateVariable(int i);
};

class SceneEdgeElasticityMarker : public SceneEdgeMarker
{
public:
    PhysicFieldBC typeX;
    PhysicFieldBC typeY;
    Value forceX;
    Value forceY;

    QString script();
    QMap<QString, QString> data();
    SceneEdgeElasticityMarker(const QString &name, PhysicFieldBC typeX, PhysicFieldBC typeY, Value forceX, Value forceY);

    int showDialog(QWidget *parent);
};

class SceneLabelElasticityMarker : public SceneLabelMarker
{
public:
    Value young_modulus;
    Value poisson_ratio;

    SceneLabelElasticityMarker(const QString &name, Value young_modulus, Value poisson_ratio);

    // Lame constant
    inline double lambda() { return (young_modulus.number * poisson_ratio.number) / ((1.0 + poisson_ratio.number) * (1.0 - 2.0*poisson_ratio.number)); }
    inline double mu() { return young_modulus.number / (2.0*(1.0 + poisson_ratio.number)); }

    QString script();
    QMap<QString, QString> data();
    int showDialog(QWidget *parent);
};

class DSceneEdgeElasticityMarker : public DSceneEdgeMarker
{
    Q_OBJECT

public:
    DSceneEdgeElasticityMarker(SceneEdgeElasticityMarker *edgeEdgeElasticityMarker, QWidget *parent);
    ~DSceneEdgeElasticityMarker();

protected:
    void createContent();

    void load();
    bool save();

private:
    QComboBox *cmbTypeX;
    QComboBox *cmbTypeY;
    SLineEditValue *txtForceX;
    SLineEditValue *txtForceY;
};

class DSceneLabelElasticityMarker : public DSceneLabelMarker
{
    Q_OBJECT

public:
    DSceneLabelElasticityMarker(QWidget *parent, SceneLabelElasticityMarker *labelElasticityMarker);
    ~DSceneLabelElasticityMarker();

protected:
    void createContent();

    void load();
    bool save();

private:
    SLineEditValue *txtYoungModulus;
    SLineEditValue *txtPoissonNumber;
};

#endif // HERMES_ELASTICITY_H
