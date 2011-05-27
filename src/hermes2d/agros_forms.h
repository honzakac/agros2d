#ifndef AGROS_FORMS_H
#define AGROS_FORMS_H


class DefaultLinearDiffusion : public WeakForm::MatrixFormVol
{
public:
    DefaultLinearDiffusion(int i, int j, scalar coeff = 1.0,
                           SymFlag sym = HERMES_SYM, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, HERMES_ANY, sym), coeff(coeff), gt(gt) { }
    DefaultLinearDiffusion(int i, int j, std::string area, scalar coeff = 1.0,
                           SymFlag sym = HERMES_SYM, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, area, sym), coeff(coeff), gt(gt) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        if (gt == HERMES_PLANAR) result = int_grad_u_grad_v<double, scalar>(n, wt, u, v);
        else {
            if (gt == HERMES_AXISYM_X) result = int_y_grad_u_grad_v<double, scalar>(n, wt, u, v, e);
            else result = int_x_grad_u_grad_v<double, scalar>(n, wt, u, v, e);
        }
        return coeff * result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result;
        if (gt == HERMES_PLANAR) result = int_grad_u_grad_v<Ord, Ord>(n, wt, u, v);
        else {
            if (gt == HERMES_AXISYM_X) result = int_y_grad_u_grad_v<Ord, Ord>(n, wt, u, v, e);
            else result = int_x_grad_u_grad_v<Ord, Ord>(n, wt, u, v, e);
        }
        return result;
    }

    // This is to make the form usable in rk_time_step().
    virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultLinearDiffusion(*this);
    }

private:
    scalar coeff;
    GeomType gt;
};

class DefaultLinearMass : public WeakForm::MatrixFormVol
{
public:
    DefaultLinearMass(int i, int j, scalar coeff = 1.0,
                      SymFlag sym = HERMES_SYM, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, HERMES_ANY, sym), coeff(coeff), gt(gt) { }
    DefaultLinearMass(int i, int j, std::string area, scalar coeff = 1.0,
                      SymFlag sym = HERMES_SYM, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, area, sym), coeff(coeff), gt(gt) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        if (gt == HERMES_PLANAR) result = int_u_v<double, scalar>(n, wt, u, v);
        else if (gt == HERMES_AXISYM_X) result = int_y_u_v<double, scalar>(n, wt, u, v, e);
        else result = int_x_u_v<double, scalar>(n, wt, u, v, e);

        return coeff * result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        if (gt == HERMES_PLANAR) result = int_u_v<Ord, Ord>(n, wt, u, v);
        else if (gt == HERMES_AXISYM_X) result = int_y_u_v<Ord, Ord>(n, wt, u, v, e);
        else result = int_x_u_v<Ord, Ord>(n, wt, u, v, e);

        return result;
    }

    // This is to make the form usable in rk_time_step().
    virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultLinearMass(*this);
    }

private:
    scalar coeff;
    GeomType gt;
};

// *******************************************************************************************************************


class DefaultLinearMagnetostatics : public WeakForm::MatrixFormVol
{
public:
    // The optional order_increase takes into account the axisymmetric part.
    DefaultLinearMagnetostatics(int i, int j, scalar coeff = 1.0,
                                SymFlag sym = HERMES_SYM, GeomType gt = HERMES_PLANAR,
                                int order_increase = 3)
        : WeakForm::MatrixFormVol(i, j, HERMES_ANY, sym), coeff(coeff), gt(gt), order_increase(order_increase) { }
    DefaultLinearMagnetostatics(int i, int j, std::string area, scalar coeff = 1.0,
                                SymFlag sym = HERMES_SYM, GeomType gt = HERMES_PLANAR, int order_increase = 3)
        : WeakForm::MatrixFormVol(i, j, area, sym), coeff(coeff), gt(gt), order_increase(order_increase) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        scalar planar_part = int_grad_u_grad_v<double, scalar>(n, wt, u, v);
        scalar axisym_part = 0;
        if (gt == HERMES_AXISYM_X)
            axisym_part = int_u_dvdy_over_y<double, scalar>(n, wt, u, v, e);
        else if (gt == HERMES_AXISYM_Y)
            axisym_part = int_u_dvdx_over_x<double, scalar>(n, wt, u, v, e);

        return coeff * (planar_part + axisym_part);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord planar_part = int_grad_u_grad_v<Ord, Ord>(n, wt, u, v);

        // This increase is for the axisymmetric part. We are not letting the
        // Ord class do it since it would automatically choose the highest order
        // due to the nonpolynomial 1/r term.
        return planar_part * Ord(order_increase);
    }

    // This is to make the form usable in rk_time_step().
    virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultLinearMagnetostatics(*this);
    }

private:
    scalar coeff;
    GeomType gt;
    int order_increase;
};

class DefaultLinearMagnetostaticsRemanence : public WeakForm::VectorFormVol
{
public:
    DefaultLinearMagnetostaticsRemanence(int i, double perm, double rem, double rem_ang, GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormVol(i), perm(perm), rem(rem), rem_ang(rem_ang), gt(gt) { }

    DefaultLinearMagnetostaticsRemanence(int i, std::string area, double perm, double rem, double rem_ang, GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormVol(i, area), perm(perm), rem(rem), rem_ang(rem_ang), gt(gt) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        for (int i = 0; i < n; i++)
            result += wt[i] * (- sin(rem_ang / 180.0 * M_PI) * v->dx[i]
                               + cos(rem_ang / 180.0 * M_PI) * v->dy[i]);

        return (gt == HERMES_PLANAR ? rem/perm : -rem/perm) * result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        for (int i = 0; i < n; i++)
            result += wt[i] * (v->dx[i] + v->dy[i]);

        return result;
    }

    // This is to make the form usable in rk_time_step().
    virtual WeakForm::VectorFormVol* clone() {
        return new DefaultLinearMagnetostaticsRemanence(*this);
    }

private:
    double perm, rem, rem_ang;
    GeomType gt;
};

class DefaultJacobianMagnetostaticsVelocity : public WeakForm::MatrixFormVol
{
public:
    DefaultJacobianMagnetostaticsVelocity(int i, int j, double gamma, double vel_x, double vel_y, double vel_ang = 0.0)
        : WeakForm::MatrixFormVol(i, j), gamma(gamma), vel_x(vel_x), vel_y(vel_y), vel_ang(vel_ang) { }

    DefaultJacobianMagnetostaticsVelocity(int i, int j, std::string area, double gamma, double vel_x, double vel_y, double vel_ang = 0.0)
        : WeakForm::MatrixFormVol(i, j, area, HERMES_NONSYM), gamma(gamma), vel_x(vel_x), vel_y(vel_y), vel_ang(vel_ang) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        for (int ii = 0; ii < n; ii++)
            result += wt[ii] * u->val[ii] * ((vel_x - e->y[ii] * vel_ang) * v->dx[ii] +
                                             (vel_y + e->x[ii] * vel_ang) * v->dy[ii]);

        return -gamma * result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        for (int ii = 0; ii < n; ii++)
            result += wt[ii] * u->val[ii] * (v->dx[ii] + v->dy[ii]);

        return result;
    }

    // This is to make the form usable in rk_time_step().
    virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultJacobianMagnetostaticsVelocity(*this);
    }

private:
    double gamma, vel_x, vel_y, vel_ang;
};

class HERMES_API DefaultResidualMagnetostaticsVelocity : public WeakForm::VectorFormVol
{
public:
    DefaultResidualMagnetostaticsVelocity(int i, std::string area,
                                          double gamma, double vel_x, double vel_y, double vel_ang = 0.0)
        : WeakForm::VectorFormVol(i, area), gamma(gamma), vel_x(vel_x), vel_y(vel_y), vel_ang(vel_ang) {}
    DefaultResidualMagnetostaticsVelocity(int i, Hermes::vector<std::string> areas,
                                          double gamma, double vel_x, double vel_y, double vel_ang = 0.0)
        : WeakForm::VectorFormVol(i, areas), gamma(gamma), vel_x(vel_x), vel_y(vel_y), vel_ang(vel_ang) {}

    ~DefaultResidualMagnetostaticsVelocity();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const
    {
        scalar result = 0;
        for (int ii = 0; ii < n; ii++)
            result += wt[ii] * u_ext[i]->val[ii] * ((vel_x - e->y[ii] * vel_ang) * v->dx[ii] +
                                                    (vel_y + e->x[ii] * vel_ang) * v->dy[ii]);
        return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const
    {
        Ord result = 0;
        for (int ii = 0; ii < n; ii++)
            result += wt[ii] * u_ext[i]->val[ii] * ((vel_x - e->y[ii] * vel_ang) * v->dx[ii] +
                                                    (vel_y + e->x[ii] * vel_ang) * v->dy[ii]);
        return result;
    }

    virtual WeakForm::VectorFormVol* clone()
    {
        return new DefaultResidualMagnetostaticsVelocity(*this);
    }

private:
    double gamma, vel_x, vel_y, vel_ang;
};

// **************************************************************************************************************

class DefaultLinearXX : public WeakForm::MatrixFormVol
{
public:
    DefaultLinearXX(unsigned int i, unsigned int j, double lambda, double mu, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu), gt(gt) { }
    DefaultLinearXX(unsigned int i, unsigned int j, std::string area, double lambda, double mu, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, area, HERMES_NONSYM), lambda(lambda), mu(mu), gt(gt) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        if (gt == HERMES_PLANAR)
            return (lambda + 2*mu) * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
                    mu * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
        else if (gt == HERMES_AXISYM_X)
        {
            Scalar result = 0;
            for (int i = 0; i < n; i++)
                result += wt[i] * (e->x[i] * lambda * (u->dx[i] * v->dx[i] +
                                                       u->val[i]/e->x[i] * v->dx[i] +
                                                       u->dx[i] * v->val[i]/e->x[i] +
                                                       1/sqr(e->x[i]) * u->val[i] * v->val[i]) +
                                   e->x[i] * mu * (2 * u->dx[i] * v->dx[i] +
                                                   2 * 1/sqr(e->x[i]) * u->val[i] * v->val[i] +
                                                   u->dy[i] * v->dy[i]));
            return result;
        }
        else
        {
            Scalar result = 0;
            for (int i = 0; i < n; i++)
                result += wt[i] * (e->x[i] * lambda * (u->dx[i] * v->dx[i] +
                                                       u->val[i]/e->x[i] * v->dx[i] +
                                                       u->dx[i] * v->val[i]/e->x[i] +
                                                       1/sqr(e->x[i]) * u->val[i] * v->val[i]) +
                                   e->x[i] * mu * (2 * u->dx[i] * v->dx[i] +
                                                   2 * 1/sqr(e->x[i]) * u->val[i] * v->val[i] +
                                                   u->dy[i] * v->dy[i]));
            return result;
        }
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // This is to make the form usable in rk_time_step().
    virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultLinearXX(*this);
    }

private:
    double lambda, mu;
    GeomType gt;
};

class DefaultLinearXY : public WeakForm::MatrixFormVol
{
public:
    DefaultLinearXY(unsigned int i, unsigned int j, double lambda, double mu, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu), gt(gt) { }
    DefaultLinearXY(unsigned int i, unsigned int j, std::string area, double lambda, double mu, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, area, HERMES_SYM), lambda(lambda), mu(mu), gt(gt) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        if (gt == HERMES_PLANAR)
            return lambda * int_dudy_dvdx<Real, Scalar>(n, wt, u, v) +
                    mu * int_dudx_dvdy<Real, Scalar>(n, wt, u, v);
        else if (gt == HERMES_AXISYM_X)
        {
            Scalar result = 0;
            for (int i = 0; i < n; i++)
                result += wt[i] * (e->y[i] * lambda * (u->dx[i] * v->dy[i] +
                                                       u->dx[i] * v->val[i]/e->y[i]) +
                                   e->y[i] * mu * u->dy[i] * v->dx[i]);
            return result;
        }
        else
        {
            Scalar result = 0;
            for (int i = 0; i < n; i++)
                result += wt[i] * (e->x[i] * lambda * (u->dy[i] * v->dx[i] +
                                                       u->dy[i] * v->val[i]/e->x[i]) +
                                   e->x[i] * mu * u->dx[i] * v->dy[i]);
            return result;
        }
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // This is to make the form usable in rk_time_step().
    virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultLinearXY(*this);
    }

private:
    double lambda, mu;
    GeomType gt;
};

class DefaultLinearYY : public WeakForm::MatrixFormVol
{
public:
    DefaultLinearYY(unsigned int i, unsigned int j, double lambda, double mu, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu), gt(gt) { }
    DefaultLinearYY(unsigned int i, unsigned int j, std::string area, double lambda, double mu, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, area, HERMES_NONSYM), lambda(lambda), mu(mu), gt(gt) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        if (gt == HERMES_PLANAR)
            return  mu * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
                    (lambda + 2*mu) * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
        else if (gt == HERMES_AXISYM_X)
        {
            Scalar result = 0;
            for (int i = 0; i < n; i++)
                result += wt[i] * (e->y[i] * lambda * (u->dx[i] * v->dx[i]) +
                                   e->y[i] * mu * (u->dy[i] * v->dy[i] +
                                                   2 * u->dx[i] * v->dx[i]));
            return result;
        }
        else
        {
            Scalar result = 0;
            for (int i = 0; i < n; i++)
                result += wt[i] * (e->x[i] * lambda * (u->dy[i] * v->dy[i]) +
                                   e->x[i] * mu * (u->dx[i] * v->dx[i] +
                                                   2 * u->dy[i] * v->dy[i]));
            return result;
        }
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // This is to make the form usable in rk_time_step().
    virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultLinearYY(*this);
    }

private:
    double lambda, mu;
    GeomType gt;
};

#endif // AGROS_FORMS_H
