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

#include "hermes_field.h"

#include "hermes_general.h"
#include "hermes_electrostatic.h"
#include "hermes_magnetic.h"
#include "hermes_heat.h"
#include "hermes_current.h"
#include "hermes_elasticity.h"
// #include "hermes_flow.h"
#include "hermes_rf.h"
#include "hermes_acoustic.h"
#include "progressdialog.h"

#include "mesh/h2d_reader.h"

class CustomInitialCondition : public ExactSolutionScalar
{
public:
    CustomInitialCondition(Mesh* mesh) : ExactSolutionScalar(mesh)
    {
    };

    virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const
    {
        dx = 0.0;
        dy = 0.0;
    };

    virtual scalar value (double x, double y) const
    {
        return 0.0;
    };

    virtual Ord ord(Ord x, Ord y) const
    {
        return Ord(1);
    }
};

double actualTime;

HermesField *hermesFieldFactory(PhysicField physicField)
{
    switch (physicField)
    {
    case PhysicField_General:
        return new HermesGeneral();
    case PhysicField_Electrostatic:
        return new HermesElectrostatic();
    case PhysicField_Magnetic:
        return new HermesMagnetic();
    case PhysicField_Heat:
        return new HermesHeat();
    case PhysicField_Current:
        return new HermesCurrent();
    case PhysicField_Elasticity:
        return new HermesElasticity();
    case PhysicField_Flow:
        // return new HermesFlow();
    case PhysicField_RF:
        return new HermesRF();
    case PhysicField_Acoustic:
        return new HermesAcoustic();
    default:
        std::cerr << "Physical field '" + QString::number(physicField).toStdString() + "' is not implemented. hermesObjectFactory()" << endl;
        throw;
        break;
    }
}

void readMeshDirtyFix()
{
    // fix precalulating matrices for mapping of curved elements

    // save locale
    char *plocale = setlocale (LC_NUMERIC, "");
    setlocale (LC_NUMERIC, "C");

    std::ostringstream os;
    os << "vertices =" << std::endl <<
          "{" << std::endl <<
          "{ 0, 0 }," << std::endl <<
          "{ 1, 0 }," << std::endl <<
          "{ 0, 1 }" << std::endl <<
          "}" << std::endl <<
          "elements =" << std::endl <<
          "{" << std::endl <<
          "{ 0, 1, 2, 0 }" << std::endl <<
          "}" << std::endl <<
          "boundaries =" << std::endl <<
          "{" << std::endl <<
          "{ 0, 1, 1 }," << std::endl <<
          "{ 1, 2, 1 }," << std::endl <<
          "{ 2, 0, 1 }" << std::endl <<
          "}" << std::endl <<
          "curves =" << std::endl <<
          "{" << std::endl <<
          "{ 1, 2, 90 }" << std::endl <<
          "}";

    Mesh mesh;
    H2DReader meshloader;
    meshloader.load_str(os.str().c_str(), &mesh);

    // set system locale
    setlocale(LC_NUMERIC, plocale);
}

Mesh *readMeshFromFile(const QString &fileName)
{
    // save locale
    char *plocale = setlocale (LC_NUMERIC, "");
    setlocale (LC_NUMERIC, "C");

    // load the mesh file
    Mesh *mesh = new Mesh();
    H2DReader meshloader;
    meshloader.load(fileName.toStdString().c_str(), mesh);

    // set system locale
    setlocale(LC_NUMERIC, plocale);

    return mesh;
}

void writeMeshFromFile(const QString &fileName, Mesh *mesh)
{
    // save locale
    char *plocale = setlocale (LC_NUMERIC, "");
    setlocale (LC_NUMERIC, "C");

    H2DReader meshloader;
    meshloader.save(fileName.toStdString().c_str(), mesh);

    // set system locale
    setlocale(LC_NUMERIC, plocale);
}

void refineMesh(Mesh *mesh, bool refineGlobal, bool refineTowardsEdge)
{
    // refine mesh - global
    if (refineGlobal)
        for (int i = 0; i < Util::scene()->problemInfo()->numberOfRefinements; i++)
            mesh->refine_all_elements(0);

    // refine mesh - boundary
    if (refineTowardsEdge)
        for (int i = 0; i < Util::scene()->edges.count(); i++)
            if (Util::scene()->edges[i]->refineTowardsEdge > 0)
                mesh->refine_towards_boundary(QString::number(((Util::scene()->edges[i]->boundary->type != PhysicFieldBC_None) ? i + 1 : -i)).toStdString(),
                                              Util::scene()->edges[i]->refineTowardsEdge);
}

// return geom type
GeomType convertProblemType(ProblemType problemType)
{
    return (problemType == ProblemType_Planar ? HERMES_PLANAR : HERMES_AXISYM_Y);
}

QList<SolutionArray *> solveSolutioArray(ProgressItemSolve *progressItemSolve,
                                         Hermes::vector<EssentialBCs> bcs,
                                         WeakFormAgros *wf)
{
    SolutionAgros solutionAgros(progressItemSolve, wf);

    QList<SolutionArray *> solutionArrayList = solutionAgros.solveSolutioArray(bcs);
    return solutionArrayList;
}

// *********************************************************************************************************************************************

SolutionAgros::SolutionAgros(ProgressItemSolve *progressItemSolve, WeakFormAgros *wf)
{
    analysisType = Util::scene()->problemInfo()->analysisType;
    polynomialOrder = Util::scene()->problemInfo()->polynomialOrder;
    adaptivityType = Util::scene()->problemInfo()->adaptivityType;
    adaptivitySteps = Util::scene()->problemInfo()->adaptivitySteps;
    adaptivityTolerance = Util::scene()->problemInfo()->adaptivityTolerance;
    adaptivityMaxDOFs = Util::scene()->problemInfo()->adaptivityMaxDOFs;
    numberOfSolution = Util::scene()->problemInfo()->hermes()->numberOfSolution();
    timeTotal = Util::scene()->problemInfo()->timeTotal.number;
    timeStep = Util::scene()->problemInfo()->timeStep.number;
    initialCondition = Util::scene()->problemInfo()->initialCondition.number;

    linearityType = Util::scene()->problemInfo()->linearityType;
    linearityNonlinearTolerance = Util::scene()->problemInfo()->linearityNonlinearTolerance;
    linearityNonlinearSteps = Util::scene()->problemInfo()->linearityNonlinearSteps;

    matrixSolver = Util::scene()->problemInfo()->matrixSolver;

    m_progressItemSolve = progressItemSolve;
    m_wf = wf;
}

QList<SolutionArray *> SolutionAgros::solveSolutioArray(Hermes::vector<EssentialBCs> bcs)
{
    QTime time;

    // solution agros array
    QList<SolutionArray *> solutionArrayList;

    // load the mesh file
    mesh = readMeshFromFile(tempProblemFileName() + ".mesh");
    refineMesh(mesh, true, true);

    // create an H1 space
    Hermes::vector<Space *> space;
    // create hermes solution array
    Hermes::vector<Solution *> solution;
    // create reference solution
    Hermes::vector<Solution *> solutionReference;

    // projection norms
    Hermes::vector<ProjNormType> projNormType;

    // prepare selector
    Hermes::vector<RefinementSelectors::Selector *> selector;

    // error marker
    isError = false;

    RefinementSelectors::Selector *select = NULL;
    switch (adaptivityType)
    {
    case AdaptivityType_H:
        select = new RefinementSelectors::HOnlySelector();
        break;
    case AdaptivityType_P:
        select = new RefinementSelectors::H1ProjBasedSelector(RefinementSelectors::H2D_P_ANISO,
                                                              Util::config()->convExp,
                                                              H2DRS_DEFAULT_ORDER);
        break;
    case AdaptivityType_HP:
        select = new RefinementSelectors::H1ProjBasedSelector(RefinementSelectors::H2D_HP_ANISO,
                                                              Util::config()->convExp,
                                                              H2DRS_DEFAULT_ORDER);
        break;
    }

    for (int i = 0; i < numberOfSolution; i++)
    {
        space.push_back(new H1Space(mesh, &bcs[i], polynomialOrder));

        // set order by element
        for (int j = 0; j < Util::scene()->labels.count(); j++)
            if (Util::scene()->labels[j]->material != Util::scene()->materials[0])
                space.at(i)->set_uniform_order(Util::scene()->labels[j]->polynomialOrder > 0 ? Util::scene()->labels[j]->polynomialOrder : polynomialOrder,
                                               QString::number(j).toStdString());

        // solution agros array
        solution.push_back(new Solution());

        if (adaptivityType != AdaptivityType_None)
        {
            // add norm
            projNormType.push_back(Util::config()->projNormType);
            // add refinement selector
            selector.push_back(select);
        }
    }

    // check for DOFs
    if (Space::get_num_dofs(space) == 0)
    {
        m_progressItemSolve->emitMessage(QObject::tr("DOF is zero"), true);
    }
    else
    {
        for (int i = 0; i < numberOfSolution; i++)
        {
            // transient
            if (analysisType == AnalysisType_Transient)
            {
                // constant initial solution
                solution.at(i)->set_const(mesh, initialCondition);
                solutionArrayList.append(solutionArray(solution.at(i)));
            }

            // nonlinear
            if ((linearityType != LinearityType_Linear) && (analysisType != AnalysisType_Transient))
            {
                solution.at(i)->set_const(mesh, 0.0);
            }
        }

        actualTime = 0.0;

        // update time function
        // Util::scene()->problemInfo()->hermes()->updateTimeFunctions(actualTime);

        // m_wf->set_current_time(actualTime);
        // m_wf->solution = solution;
        // m_wf->delete_all();
        // m_wf->registerForms();

        // emit message
        if (adaptivityType != AdaptivityType_None)
            m_progressItemSolve->emitMessage(QObject::tr("Adaptivity type: %1").arg(adaptivityTypeString(adaptivityType)), false);

        double error = 0.0;

        // solution
        /*
        int maxAdaptivitySteps = (adaptivityType == AdaptivityType_None) ? 1 : adaptivitySteps;
        int actualAdaptivitySteps = -1;
        for (int i = 0; i<maxAdaptivitySteps; i++)
        {
            // set up the solver, matrix, and rhs according to the solver selection.
            SparseMatrix *matrix = create_matrix(matrixSolver);
            Vector *rhs = create_vector(matrixSolver);
            Solver *solver = create_linear_solver(matrixSolver, matrix, rhs);

            if (adaptivityType == AdaptivityType_None)
            {
                if (analysisType != AnalysisType_Transient)
                    solve(space, solution, solver, matrix, rhs);
            }
            else
            {
                // construct globally refined reference mesh and setup reference space.
                Hermes::vector<Space *> spaceReference = *Space::construct_refined_spaces(space);

                // assemble reference problem.
                solve(spaceReference, solution, solver, matrix, rhs);

                // copy solution
                for (int j = 0; j < numberOfSolution; j++)
                {
                    solutionReference.push_back(new Solution());
                    solutionReference.at(j)->copy(solution.at(j));
                }

                if (!isError)
                {
                    // project the fine mesh solution onto the coarse mesh.
                    OGProjection::project_global(space, solutionReference, solution, matrixSolver);

                    // Calculate element errors and total error estimate.
                    Adapt adaptivity(space, projNormType);

                    // Calculate error estimate for each solution component and the total error estimate.
                    Hermes::vector<double> err_est_rel;
                    error = adaptivity.calc_err_est(solution,
                                                    solutionReference,
                                                    &err_est_rel) * 100;

                    // emit signal
                    m_progressItemSolve->emitMessage(QObject::tr("Adaptivity rel. error (step: %2/%3, DOFs: %4): %1%").
                                                     arg(error, 0, 'f', 3).
                                                     arg(i + 1).
                                                     arg(maxAdaptivitySteps).
                                                     arg(Space::get_num_dofs(space)), false, 1);
                    // add error to the list
                    m_progressItemSolve->addAdaptivityError(error, Space::get_num_dofs(space));

                    if (error < adaptivityTolerance || Space::get_num_dofs(space) >= adaptivityMaxDOFs)
                    {
                        break;
                    }
                    if (i != maxAdaptivitySteps-1) adaptivity.adapt(selector,
                                                                    Util::config()->threshold,
                                                                    Util::config()->strategy,
                                                                    Util::config()->meshRegularity);
                    actualAdaptivitySteps = i+1;
                }

                if (m_progressItemSolve->isCanceled())
                {
                    isError = true;
                    break;
                }

                // delete reference space
                for (int i = 0; i < spaceReference.size(); i++)
                {
                    delete spaceReference.at(i)->get_mesh();
                    delete spaceReference.at(i);
                }
                spaceReference.clear();

                // delete reference solution
                for (int i = 0; i < solutionReference.size(); i++)
                    delete solutionReference.at(i);
                solutionReference.clear();
            }

            // clean up.
            delete solver;
            delete matrix;
            delete rhs;
        }
        */

        int maxAdaptivitySteps = (adaptivityType == AdaptivityType_None) ? 1 : adaptivitySteps;
        int actualAdaptivitySteps = -1;

        // timesteps
        // if (!isError)
        {
            // set up the solver, matrix, and rhs according to the solver selection.
            SparseMatrix *matrix = create_matrix(matrixSolver);
            Vector *rhs = create_vector(matrixSolver);
            Solver *solver = create_linear_solver(matrixSolver, matrix, rhs);

            // allocate dp for transient solution
            // DiscreteProblem *dpTran = NULL;
            // if (analysisType == AnalysisType_Transient)
            // {
            // dpTran = new DiscreteProblem(m_wf, space);
            // }

            int timesteps = (analysisType == AnalysisType_Transient) ? floor(timeTotal/timeStep) : 1;
            for (int n = 0; n<timesteps; n++)
            {
                // set actual time
                actualTime = (n+1)*timeStep;

                // update essential bc values
                Space::update_essential_bc_values(space, actualTime);
                // update timedep values
                Util::scene()->problemInfo()->hermes()->updateTimeFunctions(actualTime);

                m_wf->set_current_time(actualTime);
                m_wf->solution = solution;
                m_wf->delete_all();
                m_wf->registerForms();

                // solve system
                solve(space, solution, solver, matrix, rhs);

                // output
                if (!isError)
                {
                    for (int i = 0; i < numberOfSolution; i++)
                        solutionArrayList.append(solutionArray(solution.at(i), space.at(i), error, actualAdaptivitySteps, (n+1)*timeStep));

                    if (analysisType == AnalysisType_Transient)
                        m_progressItemSolve->emitMessage(QObject::tr("Transient time step (%1/%2): %3 s").
                                                         arg(n+1).
                                                         arg(timesteps).
                                                         arg(actualTime, 0, 'e', 2), false, n+2);
                }
                else
                {
                    break;
                }

                if (m_progressItemSolve->isCanceled())
                {
                    isError = true;
                    break;
                }
            }

            // clean up
            if (solver) delete solver;
            if (matrix) delete matrix;
            if (rhs) delete rhs;

            // if (dpTran) delete dpTran;
        }
    }

    // delete selector
    if (select) delete select;
    selector.clear();

    // delete mesh
    delete mesh;

    // delete space
    for (unsigned int i = 0; i < space.size(); i++)
    {
        // delete space.at(i)->get_mesh();
        delete space.at(i);
    }
    space.clear();

    // delete last solution
    for (unsigned int i = 0; i < solution.size(); i++)
        delete solution.at(i);
    solution.clear();

    if (isError)
    {
        for (int i = 0; i < solutionArrayList.count(); i++)
            delete solutionArrayList.at(i);
        solutionArrayList.clear();
    }
    return solutionArrayList;
}

bool SolutionAgros::solve(Hermes::vector<Space *> space,
                          Hermes::vector<Solution *> solution,
                          Solver *solver, SparseMatrix *matrix, Vector *rhs)
{
    Hermes2D hermes2d;

    DiscreteProblem dp(m_wf, space);

    scalar *coeff_vec = new scalar[Space::get_num_dofs(space)];
    memset(coeff_vec, 0, Space::get_num_dofs(space)*sizeof(scalar));
    // CustomInitialCondition init_sln(mesh);
    // OGProjection::project_global(space, solution, &init_sln, coeff_vec);

    // Perform Newton's iteration.
    bool jacobian_changed = true;
    bool residual_as_function = false;
    double damping_coeff = 1.0;
    double max_allowed_residual_norm = 1e12;

    // Prepare solutions for measuring residual norm.
    int num_spaces = dp.get_spaces().size();
    Hermes::vector<Solution*> solutions;
    Hermes::vector<bool> dir_lift_false;
    for (int i=0; i < num_spaces; i++)
    {
        if (residual_as_function)
            solutions.push_back(new Solution());

        // No Dirichlet lifts will be considered.
        dir_lift_false.push_back(false);
    }

    // The Newton's loop.
    double residual_norm;
    int it = 1;
    while (1)
    {
        // Obtain the number of degrees of freedom.
        int ndof = dp.get_num_dofs();

        // Assemble the residual vector.
        dp.assemble(coeff_vec, NULL, rhs); // NULL = we do not want the Jacobian.

        // Calculate the l2-norm of residual vector, this is the traditional way.
        residual_norm = hermes2d.get_l2_norm(rhs);

        // Info for the user.
        if (it == 1)
            m_progressItemSolve->emitMessage(QObject::tr("Newton initial residual norm: %1").
                                             arg(residual_norm, 0, 'e', 3), false, 1);
        else
            m_progressItemSolve->emitMessage(QObject::tr("Newton iter %1, residual norm: %2").
                                             arg(it-1).
                                             arg(residual_norm, 0, 'e', 3), false, 1);

        // If maximum allowed residual norm is exceeded, fail.
        if (residual_norm > max_allowed_residual_norm)
        {
            m_progressItemSolve->emitMessage(QObject::tr("Current residual norm: %1").
                                             arg(residual_norm, 0, 'e', 3), false, 1);
            m_progressItemSolve->emitMessage(QObject::tr("Maximum allowed residual norm: %1").
                                             arg(max_allowed_residual_norm, 0, 'e', 3), false, 1);
            m_progressItemSolve->emitMessage(QObject::tr("Newton solve not successful."), true, 1);

            isError = true;
            break;
        }

        // If residual norm is within tolerance, or the maximum number
        // of iteration has been reached, then quit.
        if ((residual_norm < Util::scene()->problemInfo()->linearityNonlinearTolerance || it > Util::scene()->problemInfo()->linearityNonlinearSteps) && it > 1)
            break;

        // If Jacobian changed, assemble the matrix.
        if (jacobian_changed)
            dp.assemble(coeff_vec, matrix, NULL); // NULL = we do not want the rhs.

        // Multiply the residual vector with -1 since the matrix
        // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
        rhs->change_sign();

        // Solve the linear system.
        if(!solver->solve())
        {
            m_progressItemSolve->emitMessage(QObject::tr("Matrix solver failed."), true, 1);
            isError = true;
            break;
        }

        // Add \deltaY^{n+1} to Y^n.
        for (int i = 0; i < ndof; i++)
            coeff_vec[i] += damping_coeff * solver->get_solution()[i];

        it++;
    }

    if (it >= Util::scene()->problemInfo()->linearityNonlinearSteps)
    {
        m_progressItemSolve->emitMessage(QObject::tr("Maximum allowed number of Newton iterations exceeded."), true, 1);
        isError = true;
    }

    // Translate the resulting coefficient vector into the Solution sln.
    if (!isError)
        Solution::vector_to_solutions(coeff_vec, space, solution);

    // Clean up.
    delete [] coeff_vec;

    return !isError;
}

SolutionArray *SolutionAgros::solutionArray(Solution *sln, Space *space, double adaptiveError, double adaptiveSteps, double time)
{
    SolutionArray *solution = new SolutionArray();
    solution->order = new Orderizer();
    if (space) solution->order->process_space(space);
    solution->sln = new Solution();
    if (sln) solution->sln->copy(sln);
    solution->adaptiveError = adaptiveError;
    solution->adaptiveSteps = adaptiveSteps;
    solution->time = time;

    return solution;
}

// *********************************************************************************************************************************************

ViewScalarFilter::ViewScalarFilter(Hermes::vector<MeshFunction *> sln, PhysicFieldVariable physicFieldVariable, PhysicFieldVariableComp physicFieldVariableComp)
    : Filter(sln)
{
    m_physicFieldVariable = physicFieldVariable;
    m_physicFieldVariableComp = physicFieldVariableComp;
}

double ViewScalarFilter::get_pt_value(double x, double y, int item)
{
    return 0.0;
}

void ViewScalarFilter::precalculate(int order, int mask)
{
    Quad2D* quad = quads[cur_quad];
    int np = quad->get_num_points(order);
    node = new_node(H2D_FN_DEFAULT, np);

    if (sln[0])
    {
        sln[0]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
        sln[0]->get_dx_dy_values(dudx1, dudy1);
        value1 = sln[0]->get_fn_values();
    }

    if (num >= 2 && sln[1])
    {
        sln[1]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
        sln[1]->get_dx_dy_values(dudx2, dudy2);
        value2 = sln[1]->get_fn_values();
    }

    if (num >= 3 && sln[2])
    {
        sln[2]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
        sln[2]->get_dx_dy_values(dudx3, dudy3);
        value3 = sln[2]->get_fn_values();
    }

    update_refmap();

    x = refmap->get_phys_x(order);
    y = refmap->get_phys_y(order);
    Element *e = refmap->get_active_element();

    material = Util::scene()->labels[Util::scene()->sceneSolution()->agrosMaterialMarker(e->marker)]->material;

    for (int i = 0; i < np; i++)
    {
        calculateVariable(i);
    }

    if (nodes->present(order))
    {
        assert(nodes->get(order) == cur_node);
        ::free(nodes->get(order));
    }
    nodes->add(node, order);
    cur_node = node;
}
