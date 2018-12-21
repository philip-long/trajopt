#ifndef KINEMATICS_TERMS_HPP
#define KINEMATICS_TERMS_HPP

#pragma once

#include "sco/modeling.hpp"
#include "sco/modeling_utils.hpp"
#include "sco/sco_fwd.hpp"
#include <Eigen/Core>
#include "trajopt/common.hpp"
#include <openrave/openrave.h>
namespace trajopt {

using namespace sco;
typedef BasicArray<Var> VarArray;

#if 0
void makeTrajVariablesAndBounds(int n_steps, const RobotAndDOF& manip, OptProb& prob_out, VarArray& vars_out);

class FKFunc {
public:
  virtual OpenRAVE::Transform operator()(const VectorXd& x) const = 0;
  virtual ~FKFunc() {}
};

class FKPositionJacobian {
public:
  virtual Eigen::MatrixXd operator()(const VectorXd& x) const = 0;
  virtual ~FKPositionJacobian() {}
};
#endif


struct CartPoseErrCalculator : public VectorOfVector {
  OR::Transform pose_inv_;
  ConfigurationPtr manip_;
  OR::KinBody::LinkPtr link_;
  Vector3d offset_;
  CartPoseErrCalculator(const OR::Transform& pose, ConfigurationPtr manip, OR::KinBody::LinkPtr link,Eigen::Vector3d offset) :
    pose_inv_(pose.inverse()),
    manip_(manip),
    link_(link),
    offset_(offset){}
  VectorXd operator()(const VectorXd& dof_vals) const;
};

struct CartPoseErrorPlotter : public Plotter {
  boost::shared_ptr<void> m_calc; //actually points to a CartPoseErrCalculator = CartPoseCost::f_
  VarVector m_vars;
  CartPoseErrorPlotter(boost::shared_ptr<void> calc, const VarVector& vars) : m_calc(calc), m_vars(vars) {}
  void Plot(const DblVec& x, OR::EnvironmentBase& env, std::vector<OR::GraphHandlePtr>& handles);
};

struct CartDDPoseErrCalculator : public VectorOfVector {
  ConfigurationPtr manip_;
  OR::KinBody::LinkPtr link_;
  Vector3d offset_;
  CartDDPoseErrCalculator(ConfigurationPtr manip, OR::KinBody::LinkPtr link,Eigen::Vector3d offset) :
  manip_(manip),
  link_(link),
  offset_(offset){}
  VectorXd operator()(const VectorXd& dof_vals) const;
};

// Philip Adding this pinv function, it looks wierd because trajopt compliation
// is not compatible with C++ 11,
// the original is herehttps://gist.github.com/javidcf/25066cf85e71105d57b6
//--------------------KINEMATIC FUNCTIONS FOR VALKYRIE--------------------------//
// This needs to be moved somewhere else---------------------------//
template <class MatT>
Eigen::Matrix<typename MatT::Scalar, MatT::ColsAtCompileTime, MatT::RowsAtCompileTime>
pseudoinverse(const MatT &mat, double tolerance = 0.0001) // choose appropriately
{
    typedef typename MatT::Scalar Scalar;

    Eigen::JacobiSVD<Eigen::Matrix<Scalar,MatT::RowsAtCompileTime,MatT::ColsAtCompileTime> > svd = mat.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix<Scalar,MatT::RowsAtCompileTime,1>  singularValues = svd.singularValues();
    Eigen::Matrix<Scalar, MatT::ColsAtCompileTime, MatT::RowsAtCompileTime> singularValuesInv(mat.cols(), mat.rows());
    singularValuesInv.setZero();
    for (unsigned int i = 0; i < singularValues.size(); ++i) {
        if (singularValues(i) > tolerance)
        {
            singularValuesInv(i, i) = 1.0/ singularValues(i);
        }
        else
        {
            singularValuesInv(i, i) =0.0;
        }
    }
    return svd.matrixV() * singularValuesInv * svd.matrixU().adjoint();
}

void getJacobianTorsoRightpalm(std::vector<double> q,Eigen::Matrix<double,6,7> &J);
void getJacobianTorsoRightelbowpitchlink(std::vector<double> q,Eigen::Matrix<double,6,7> &J);

struct VirtualGuideErrCalculator : public VectorOfVector {
  ConfigurationPtr manip_;
  OR::EnvironmentBasePtr env_;
  OR::KinBody::LinkPtr link_;
  Vector3d translational_deviation_,angular_deviation_;


  VirtualGuideErrCalculator(OR::EnvironmentBasePtr env,
                            ConfigurationPtr manip,
                            OR::KinBody::LinkPtr link,
                            Eigen::Vector3d translational_deviation,
                            Eigen::Vector3d angular_deviation) :
  env_(env),
  manip_(manip),
  link_(link),
  translational_deviation_(translational_deviation),
  angular_deviation_(angular_deviation){}


  // I think I need to add the object name.
  VectorXd operator()(const VectorXd& dof_vals) const;
};

struct CartPoseConstraintCalculator : public VectorOfVector {
	OR::Vector plane1_;
	OR::Vector plane2_;
	OR::Vector normal_;
	ConfigurationPtr manip_;
	OR::KinBody::LinkPtr link_;
	CartPoseConstraintCalculator(const OR::Vector plane1, const OR::Vector plane2, ConfigurationPtr manip, OR::KinBody::LinkPtr link) :
	plane1_(plane1),
	plane2_(plane2),
    manip_(manip),
	link_(link)
	{normal_ = (plane1 - plane2).normalize();}
	VectorXd operator()(const VectorXd& dof_vals) const;

};

struct CartVelJacCalculator : MatrixOfVector {
  ConfigurationPtr manip_;
  KinBody::LinkPtr link_;
  double limit_;
  CartVelJacCalculator(ConfigurationPtr manip, KinBody::LinkPtr link, double limit) :
    manip_(manip), link_(link), limit_(limit) {}

  MatrixXd operator()(const VectorXd& dof_vals) const;
};

struct CartVelCalculator : VectorOfVector {
  ConfigurationPtr manip_;
  KinBody::LinkPtr link_;
  double limit_;
  CartVelCalculator(ConfigurationPtr manip, KinBody::LinkPtr link, double limit) :
    manip_(manip), link_(link), limit_(limit) {}

  VectorXd operator()(const VectorXd& dof_vals) const;
};

#if 0
class CartPoseCost : public CostFromErrFunc {
public:
  CartPoseCost(const VarVector& vars, const OR::Transform& pose, RobotAndDOFPtr manip, KinBody::LinkPtr link, const VectorXd& coeffs);
};

class CartPoseConstraint : public ConstraintFromFunc {
public:
  CartPoseConstraint(const VarVector& vars, const OR::Transform& pose, RobotAndDOFPtr manip, KinBody::LinkPtr link, const VectorXd& coeffs);
};

class CartVelConstraint : public ConstraintFromFunc {
public:
  CartVelConstraint(const VarVector& step0vars, const VarVector& step1vars, RobotAndDOFPtr manip, KinBody::LinkPtr link, double distlimit);
};
#endif



}

#endif // KINEMATIC_TERMS_HPP
