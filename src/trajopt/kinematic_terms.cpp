#include "sco/expr_ops.hpp"
#include "sco/modeling_utils.hpp"
#include "trajopt/kinematic_terms.hpp"
#include "trajopt/rave_utils.hpp"
#include "trajopt/utils.hpp"
#include "utils/eigen_conversions.hpp"
#include "utils/eigen_slicing.hpp"
#include "utils/logging.hpp"
#include "utils/stl_to_string.hpp"
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <Eigen/Geometry>
#include <iostream>

using namespace std;
using namespace sco;
using namespace Eigen;
using namespace util;


namespace {
  
#if 0
Vector3d rotVec(const Matrix3d& m) {
  Quaterniond q; q = m;
  return Vector3d(q.x(), q.y(), q.z());
}
#endif
inline Vector3d rotVec(const OpenRAVE::Vector& q) {
  return Vector3d(q[1], q[2], q[3]);
}

#if 0
VectorXd concat(const VectorXd& a, const VectorXd& b) {
  VectorXd out(a.size()+b.size());
  out.topRows(a.size()) = a;
  out.middleRows(a.size(), b.size()) = b;
  return out;
}

template <typename T>
vector<T> concat(const vector<T>& a, const vector<T>& b) {
  vector<T> out;
  vector<int> x;
  out.insert(out.end(), a.begin(), a.end());
  out.insert(out.end(), b.begin(), b.end());
  return out;
}
#endif

}

namespace trajopt {




// CostPtr ConstructCost(VectorOfVectorPtr err_calc, const VarVector& vars, const VectorXd& coeffs, PenaltyType type, const string& name) {
//   return CostPtr(new CostFromErrFunc(err_calc), vars, coeffs, type, name);
// }
  
  
VectorXd CartPoseErrCalculator::operator()(const VectorXd& dof_vals) const {
  manip_->SetDOFValues(toDblVec(dof_vals));
  OR::Transform newpose = link_->GetTransform();
  newpose.trans = link_->GetTransform()*OR::Vector(offset_[0],offset_[1],offset_[2],0);

  OR::Transform pose_err = pose_inv_ * newpose;
  VectorXd err = concat(rotVec(pose_err.rot), toVector3d(pose_err.trans));
  return err;  
}

VectorXd CartDDPoseErrCalculator::operator()(const VectorXd& dof_vals) const {
  int n_dof = manip_->GetDOF();
  manip_->SetDOFValues(toDblVec(dof_vals.topRows(n_dof)));
  OR::Transform pose = link_->GetTransform();
  OR::Transform pose_inv_ = pose.inverse();
  manip_->SetDOFValues(toDblVec(dof_vals.bottomRows(n_dof)));
  OR::Transform newpose = link_->GetTransform();
  newpose.trans = link_->GetTransform()*OR::Vector(offset_[0],offset_[1],offset_[2],0);

  OR::Transform pose_err = pose_inv_ * newpose;
  VectorXd err = concat(rotVec(pose_err.rot), toVector3d(pose_err.trans));
  return err;  
}

// Philip addingsdvsdvds


VectorXd VirtualGuideErrCalculator::operator()(const VectorXd& dof_vals) const {

  int n_dof = manip_->GetDOF();
  double k(0.0);
  std::cout<<"n_dof"<<n_dof<<std::endl;
  manip_->SetDOFValues(toDblVec(dof_vals.topRows(n_dof)));
  DblVec joint_vals=manip_->GetDOFValues();
  // The first element of the right arm joint is 26
  // last element is 32
  DblVec right_joint_val(&joint_vals[25],&joint_vals[32]);
  Eigen::Matrix<double,6,7> b_J_e,b_J_ep;
  Eigen::Matrix<double,7,6> b_J_e_inv;
  getJacobianTorsoRightpalm(right_joint_val,b_J_e);
  getJacobianTorsoRightelbowpitchlink(right_joint_val,b_J_ep);
  b_J_e_inv=pseudoinverse(b_J_e);

  // I need to find the closest point distance and location
  // Project Jacobian along the normal to this point


  OR::Transform pose = link_->GetTransform();
  OR::Transform pose_inv_ = pose.inverse();
  manip_->SetDOFValues(toDblVec(dof_vals.bottomRows(n_dof)));
  OR::Transform newpose = link_->GetTransform();
  OR::Transform pose_err = pose_inv_ * newpose;
  VectorXd err = concat(rotVec(pose_err.rot), toVector3d(pose_err.trans));
  return err;
}




void getJacobianTorsoRightpalm(std::vector<double> q,Eigen::Matrix<double,6,7> &J)
 {
 J.setZero();
double rightShoulderPitch=q[0];
double rightShoulderRoll=q[1];
double rightShoulderYaw=q[2];
double rightElbowPitch=q[3];
double rightForearmYaw=q[4];
double rightWristRoll=q[5];
double rightWristPitch=q[6];


J(0,0)=-0.2871*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*sin(rightElbowPitch) + 0.0254*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*cos(rightElbowPitch) - 0.0254*sin(rightElbowPitch)*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 0.0254*sin(rightShoulderPitch)*cos(rightShoulderYaw) - 0.2871*sin(rightShoulderRoll)*cos(rightElbowPitch)*cos(rightShoulderPitch) - 0.33*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 0.0254*sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll);
J(0,1)=1.0*(0.2871*sin(rightElbowPitch)*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.0254*sin(rightElbowPitch)*cos(rightShoulderRoll) - 0.0254*sin(rightShoulderRoll)*sin(rightShoulderYaw)*cos(rightElbowPitch) + 0.0254*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.2871*cos(rightElbowPitch)*cos(rightShoulderRoll) - 0.33*cos(rightShoulderRoll))*sin(rightShoulderPitch);
J(0,2)=-1.0*(0.2871*sin(rightElbowPitch)*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.0254*sin(rightElbowPitch)*cos(rightShoulderRoll) - 0.0254*sin(rightShoulderRoll)*sin(rightShoulderYaw)*cos(rightElbowPitch) + 0.0254*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.2871*cos(rightElbowPitch)*cos(rightShoulderRoll) - 0.33*cos(rightShoulderRoll))*sin(rightShoulderRoll)*cos(rightShoulderPitch) + 1.0*(-0.2871*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*sin(rightElbowPitch) + 0.0254*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*cos(rightElbowPitch) - 0.0254*sin(rightElbowPitch)*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 0.0254*sin(rightShoulderPitch)*cos(rightShoulderYaw) - 0.2871*sin(rightShoulderRoll)*cos(rightElbowPitch)*cos(rightShoulderPitch) - 0.33*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 0.0254*sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*cos(rightShoulderRoll);
J(0,3)=(1.0*sin(rightShoulderPitch)*sin(rightShoulderYaw) - 1.0*cos(rightShoulderPitch)*cos(rightShoulderRoll)*cos(rightShoulderYaw))*(0.2871*sin(rightElbowPitch)*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.0254*sin(rightElbowPitch)*cos(rightShoulderRoll) - 0.0254*sin(rightShoulderRoll)*sin(rightShoulderYaw)*cos(rightElbowPitch) - 0.2871*cos(rightElbowPitch)*cos(rightShoulderRoll)) - 1.0*(-0.2871*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*sin(rightElbowPitch) + 0.0254*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*cos(rightElbowPitch) - 0.0254*sin(rightElbowPitch)*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 0.2871*sin(rightShoulderRoll)*cos(rightElbowPitch)*cos(rightShoulderPitch))*sin(rightShoulderRoll)*cos(rightShoulderYaw);
J(0,4)=(-1.0*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*sin(rightElbowPitch) - 1.0*sin(rightShoulderRoll)*cos(rightElbowPitch)*cos(rightShoulderPitch))*(0.2871*sin(rightElbowPitch)*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.2871*cos(rightElbowPitch)*cos(rightShoulderRoll)) + (-0.2871*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*sin(rightElbowPitch) - 0.2871*sin(rightShoulderRoll)*cos(rightElbowPitch)*cos(rightShoulderPitch))*(-1.0*sin(rightElbowPitch)*sin(rightShoulderRoll)*sin(rightShoulderYaw) + 1.0*cos(rightElbowPitch)*cos(rightShoulderRoll));
J(0,5)=0;
J(0,6)=0;
J(1,0)=0;
J(1,1)=-1.0*(-0.2871*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*sin(rightElbowPitch) + 0.0254*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*cos(rightElbowPitch) - 0.0254*sin(rightElbowPitch)*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 0.0254*sin(rightShoulderPitch)*cos(rightShoulderYaw) - 0.2871*sin(rightShoulderRoll)*cos(rightElbowPitch)*cos(rightShoulderPitch) - 0.33*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 0.0254*sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*cos(rightShoulderPitch) - 1.0*(-0.2871*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightElbowPitch) + 0.0254*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*cos(rightElbowPitch) - 0.0254*sin(rightElbowPitch)*sin(rightShoulderPitch)*sin(rightShoulderRoll) - 0.2871*sin(rightShoulderPitch)*sin(rightShoulderRoll)*cos(rightElbowPitch) - 0.33*sin(rightShoulderPitch)*sin(rightShoulderRoll) - 0.0254*sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) + 0.0254*cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightShoulderPitch);
J(1,2)=-1.0*(-0.2871*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*sin(rightElbowPitch) + 0.0254*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*cos(rightElbowPitch) - 0.0254*sin(rightElbowPitch)*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 0.0254*sin(rightShoulderPitch)*cos(rightShoulderYaw) - 0.2871*sin(rightShoulderRoll)*cos(rightElbowPitch)*cos(rightShoulderPitch) - 0.33*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 0.0254*sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*sin(rightShoulderPitch)*sin(rightShoulderRoll) + 1.0*(-0.2871*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightElbowPitch) + 0.0254*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*cos(rightElbowPitch) - 0.0254*sin(rightElbowPitch)*sin(rightShoulderPitch)*sin(rightShoulderRoll) - 0.2871*sin(rightShoulderPitch)*sin(rightShoulderRoll)*cos(rightElbowPitch) - 0.33*sin(rightShoulderPitch)*sin(rightShoulderRoll) - 0.0254*sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) + 0.0254*cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightShoulderRoll)*cos(rightShoulderPitch);
J(1,3)=(-1.0*sin(rightShoulderPitch)*sin(rightShoulderYaw) + 1.0*cos(rightShoulderPitch)*cos(rightShoulderRoll)*cos(rightShoulderYaw))*(-0.2871*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightElbowPitch) + 0.0254*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*cos(rightElbowPitch) - 0.0254*sin(rightElbowPitch)*sin(rightShoulderPitch)*sin(rightShoulderRoll) - 0.2871*sin(rightShoulderPitch)*sin(rightShoulderRoll)*cos(rightElbowPitch)) + (-1.0*sin(rightShoulderPitch)*cos(rightShoulderRoll)*cos(rightShoulderYaw) - 1.0*sin(rightShoulderYaw)*cos(rightShoulderPitch))*(-0.2871*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*sin(rightElbowPitch) + 0.0254*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*cos(rightElbowPitch) - 0.0254*sin(rightElbowPitch)*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 0.2871*sin(rightShoulderRoll)*cos(rightElbowPitch)*cos(rightShoulderPitch));
J(1,4)=(-0.2871*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*sin(rightElbowPitch) - 0.2871*sin(rightShoulderRoll)*cos(rightElbowPitch)*cos(rightShoulderPitch))*(-1.0*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightElbowPitch) - 1.0*sin(rightShoulderPitch)*sin(rightShoulderRoll)*cos(rightElbowPitch)) + (1.0*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*sin(rightElbowPitch) + 1.0*sin(rightShoulderRoll)*cos(rightElbowPitch)*cos(rightShoulderPitch))*(-0.2871*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightElbowPitch) - 0.2871*sin(rightShoulderPitch)*sin(rightShoulderRoll)*cos(rightElbowPitch));
J(1,5)=0;
J(1,6)=0;
J(2,0)=0.2871*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightElbowPitch) - 0.0254*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*cos(rightElbowPitch) + 0.0254*sin(rightElbowPitch)*sin(rightShoulderPitch)*sin(rightShoulderRoll) + 0.2871*sin(rightShoulderPitch)*sin(rightShoulderRoll)*cos(rightElbowPitch) + 0.33*sin(rightShoulderPitch)*sin(rightShoulderRoll) + 0.0254*sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - 0.0254*cos(rightShoulderPitch)*cos(rightShoulderYaw);
J(2,1)=1.0*(0.2871*sin(rightElbowPitch)*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.0254*sin(rightElbowPitch)*cos(rightShoulderRoll) - 0.0254*sin(rightShoulderRoll)*sin(rightShoulderYaw)*cos(rightElbowPitch) + 0.0254*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.2871*cos(rightElbowPitch)*cos(rightShoulderRoll) - 0.33*cos(rightShoulderRoll))*cos(rightShoulderPitch);
J(2,2)=1.0*(0.2871*sin(rightElbowPitch)*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.0254*sin(rightElbowPitch)*cos(rightShoulderRoll) - 0.0254*sin(rightShoulderRoll)*sin(rightShoulderYaw)*cos(rightElbowPitch) + 0.0254*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.2871*cos(rightElbowPitch)*cos(rightShoulderRoll) - 0.33*cos(rightShoulderRoll))*sin(rightShoulderPitch)*sin(rightShoulderRoll) - 1.0*(-0.2871*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightElbowPitch) + 0.0254*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*cos(rightElbowPitch) - 0.0254*sin(rightElbowPitch)*sin(rightShoulderPitch)*sin(rightShoulderRoll) - 0.2871*sin(rightShoulderPitch)*sin(rightShoulderRoll)*cos(rightElbowPitch) - 0.33*sin(rightShoulderPitch)*sin(rightShoulderRoll) - 0.0254*sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) + 0.0254*cos(rightShoulderPitch)*cos(rightShoulderYaw))*cos(rightShoulderRoll);
J(2,3)=(1.0*sin(rightShoulderPitch)*cos(rightShoulderRoll)*cos(rightShoulderYaw) + 1.0*sin(rightShoulderYaw)*cos(rightShoulderPitch))*(0.2871*sin(rightElbowPitch)*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.0254*sin(rightElbowPitch)*cos(rightShoulderRoll) - 0.0254*sin(rightShoulderRoll)*sin(rightShoulderYaw)*cos(rightElbowPitch) - 0.2871*cos(rightElbowPitch)*cos(rightShoulderRoll)) + 1.0*(-0.2871*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightElbowPitch) + 0.0254*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*cos(rightElbowPitch) - 0.0254*sin(rightElbowPitch)*sin(rightShoulderPitch)*sin(rightShoulderRoll) - 0.2871*sin(rightShoulderPitch)*sin(rightShoulderRoll)*cos(rightElbowPitch))*sin(rightShoulderRoll)*cos(rightShoulderYaw);
J(2,4)=(-0.2871*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightElbowPitch) - 0.2871*sin(rightShoulderPitch)*sin(rightShoulderRoll)*cos(rightElbowPitch))*(1.0*sin(rightElbowPitch)*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 1.0*cos(rightElbowPitch)*cos(rightShoulderRoll)) + (1.0*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightElbowPitch) + 1.0*sin(rightShoulderPitch)*sin(rightShoulderRoll)*cos(rightElbowPitch))*(0.2871*sin(rightElbowPitch)*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.2871*cos(rightElbowPitch)*cos(rightShoulderRoll));
J(2,5)=0;
J(2,6)=0;
J(3,0)=0;
J(3,1)=1.0*cos(rightShoulderPitch);
J(3,2)=1.0*sin(rightShoulderPitch)*sin(rightShoulderRoll);
J(3,3)=1.0*sin(rightShoulderPitch)*cos(rightShoulderRoll)*cos(rightShoulderYaw) + 1.0*sin(rightShoulderYaw)*cos(rightShoulderPitch);
J(3,4)=1.0*(sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightElbowPitch) + 1.0*sin(rightShoulderPitch)*sin(rightShoulderRoll)*cos(rightElbowPitch);
J(3,5)=-1.0*((sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*cos(rightElbowPitch) - sin(rightElbowPitch)*sin(rightShoulderPitch)*sin(rightShoulderRoll))*cos(rightForearmYaw) - 1.0*(sin(rightShoulderPitch)*cos(rightShoulderRoll)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch))*sin(rightForearmYaw);
J(3,6)=-1.0*(((sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*cos(rightElbowPitch) - sin(rightElbowPitch)*sin(rightShoulderPitch)*sin(rightShoulderRoll))*sin(rightForearmYaw) - (sin(rightShoulderPitch)*cos(rightShoulderRoll)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch))*cos(rightForearmYaw))*cos(rightWristRoll) - 1.0*((sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightElbowPitch) + sin(rightShoulderPitch)*sin(rightShoulderRoll)*cos(rightElbowPitch))*sin(rightWristRoll);
J(4,0)=1.00000000000000;
J(4,1)=0;
J(4,2)=1.0*cos(rightShoulderRoll);
J(4,3)=-1.0*sin(rightShoulderRoll)*cos(rightShoulderYaw);
J(4,4)=-1.0*sin(rightElbowPitch)*sin(rightShoulderRoll)*sin(rightShoulderYaw) + 1.0*cos(rightElbowPitch)*cos(rightShoulderRoll);
J(4,5)=1.0*(sin(rightElbowPitch)*cos(rightShoulderRoll) + sin(rightShoulderRoll)*sin(rightShoulderYaw)*cos(rightElbowPitch))*cos(rightForearmYaw) + 1.0*sin(rightForearmYaw)*sin(rightShoulderRoll)*cos(rightShoulderYaw);
J(4,6)=1.0*((sin(rightElbowPitch)*cos(rightShoulderRoll) + sin(rightShoulderRoll)*sin(rightShoulderYaw)*cos(rightElbowPitch))*sin(rightForearmYaw) - sin(rightShoulderRoll)*cos(rightForearmYaw)*cos(rightShoulderYaw))*cos(rightWristRoll) + 1.0*(sin(rightElbowPitch)*sin(rightShoulderRoll)*sin(rightShoulderYaw) - cos(rightElbowPitch)*cos(rightShoulderRoll))*sin(rightWristRoll);
J(5,0)=0;
J(5,1)=-1.0*sin(rightShoulderPitch);
J(5,2)=1.0*sin(rightShoulderRoll)*cos(rightShoulderPitch);
J(5,3)=-1.0*sin(rightShoulderPitch)*sin(rightShoulderYaw) + 1.0*cos(rightShoulderPitch)*cos(rightShoulderRoll)*cos(rightShoulderYaw);
J(5,4)=1.0*(sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*sin(rightElbowPitch) + 1.0*sin(rightShoulderRoll)*cos(rightElbowPitch)*cos(rightShoulderPitch);
J(5,5)=-1.0*((sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*cos(rightElbowPitch) - sin(rightElbowPitch)*sin(rightShoulderRoll)*cos(rightShoulderPitch))*cos(rightForearmYaw) + 1.0*(sin(rightShoulderPitch)*sin(rightShoulderYaw) - cos(rightShoulderPitch)*cos(rightShoulderRoll)*cos(rightShoulderYaw))*sin(rightForearmYaw);
J(5,6)=-1.0*(((sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*cos(rightElbowPitch) - sin(rightElbowPitch)*sin(rightShoulderRoll)*cos(rightShoulderPitch))*sin(rightForearmYaw) + (sin(rightShoulderPitch)*sin(rightShoulderYaw) - cos(rightShoulderPitch)*cos(rightShoulderRoll)*cos(rightShoulderYaw))*cos(rightForearmYaw))*cos(rightWristRoll) - 1.0*((sin(rightShoulderPitch)*cos(rightShoulderYaw) + sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*sin(rightElbowPitch) + sin(rightShoulderRoll)*cos(rightElbowPitch)*cos(rightShoulderPitch))*sin(rightWristRoll);
}


void getJacobianTorsoRightelbowpitchlink(std::vector<double> q,Eigen::Matrix<double,6,7> &J)
 {
 J.setZero();
double rightShoulderPitch=q[0];
double rightShoulderRoll=q[1];
double rightShoulderYaw=q[2];
double rightElbowPitch=q[3];
double rightForearmYaw=q[4];
double rightWristRoll=q[5];
double rightWristPitch=q[6];


J(0,0)=-0.0254*sin(rightShoulderPitch)*cos(rightShoulderYaw) - 0.33*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 0.0254*sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll);
J(0,1)=1.0*(0.0254*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.33*cos(rightShoulderRoll))*sin(rightShoulderPitch);
J(0,2)=-1.0*(0.0254*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.33*cos(rightShoulderRoll))*sin(rightShoulderRoll)*cos(rightShoulderPitch) + 1.0*(-0.0254*sin(rightShoulderPitch)*cos(rightShoulderYaw) - 0.33*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 0.0254*sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*cos(rightShoulderRoll);
J(0,3)=0;
J(1,0)=0;
J(1,1)=-1.0*(-0.33*sin(rightShoulderPitch)*sin(rightShoulderRoll) - 0.0254*sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) + 0.0254*cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightShoulderPitch) - 1.0*(-0.0254*sin(rightShoulderPitch)*cos(rightShoulderYaw) - 0.33*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 0.0254*sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*cos(rightShoulderPitch);
J(1,2)=1.0*(-0.33*sin(rightShoulderPitch)*sin(rightShoulderRoll) - 0.0254*sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) + 0.0254*cos(rightShoulderPitch)*cos(rightShoulderYaw))*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 1.0*(-0.0254*sin(rightShoulderPitch)*cos(rightShoulderYaw) - 0.33*sin(rightShoulderRoll)*cos(rightShoulderPitch) - 0.0254*sin(rightShoulderYaw)*cos(rightShoulderPitch)*cos(rightShoulderRoll))*sin(rightShoulderPitch)*sin(rightShoulderRoll);
J(1,3)=0;
J(2,0)=0.33*sin(rightShoulderPitch)*sin(rightShoulderRoll) + 0.0254*sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) - 0.0254*cos(rightShoulderPitch)*cos(rightShoulderYaw);
J(2,1)=1.0*(0.0254*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.33*cos(rightShoulderRoll))*cos(rightShoulderPitch);
J(2,2)=1.0*(0.0254*sin(rightShoulderRoll)*sin(rightShoulderYaw) - 0.33*cos(rightShoulderRoll))*sin(rightShoulderPitch)*sin(rightShoulderRoll) - 1.0*(-0.33*sin(rightShoulderPitch)*sin(rightShoulderRoll) - 0.0254*sin(rightShoulderPitch)*sin(rightShoulderYaw)*cos(rightShoulderRoll) + 0.0254*cos(rightShoulderPitch)*cos(rightShoulderYaw))*cos(rightShoulderRoll);
J(2,3)=0;
J(3,0)=0;
J(3,1)=1.0*cos(rightShoulderPitch);
J(3,2)=1.0*sin(rightShoulderPitch)*sin(rightShoulderRoll);
J(3,3)=1.0*sin(rightShoulderPitch)*cos(rightShoulderRoll)*cos(rightShoulderYaw) + 1.0*sin(rightShoulderYaw)*cos(rightShoulderPitch);
J(4,0)=1.00000000000000;
J(4,1)=0;
J(4,2)=1.0*cos(rightShoulderRoll);
J(4,3)=-1.0*sin(rightShoulderRoll)*cos(rightShoulderYaw);
J(5,0)=0;
J(5,1)=-1.0*sin(rightShoulderPitch);
J(5,2)=1.0*sin(rightShoulderRoll)*cos(rightShoulderPitch);
J(5,3)=-1.0*sin(rightShoulderPitch)*sin(rightShoulderYaw) + 1.0*cos(rightShoulderPitch)*cos(rightShoulderRoll)*cos(rightShoulderYaw);
}




#if 0
CartPoseCost::CartPoseCost(const VarVector& vars, const OR::Transform& pose, RobotAndDOFPtr manip, KinBody::LinkPtr link, const VectorXd& coeffs) :
    CostFromErrFunc(VectorOfVectorPtr(new CartPoseErrCalculator(pose, manip, link)), vars, coeffs, ABS, "CartPose")
{}
CartPoseConstraint::CartPoseConstraint(const VarVector& vars, const OR::Transform& pose,
    RobotAndDOFPtr manip, KinBody::LinkPtr link, const VectorXd& coeffs) :
    ConstraintFromFunc(VectorOfVectorPtr(new CartPoseErrCalculator(pose, manip, link)), vars, coeffs, EQ, "CartPose")
{}
#endif

void CartPoseErrorPlotter::Plot(const DblVec& x, OR::EnvironmentBase& env, std::vector<OR::GraphHandlePtr>& handles) {
  CartPoseErrCalculator* calc = static_cast<CartPoseErrCalculator*>(m_calc.get());
  DblVec dof_vals = getDblVec(x, m_vars);
  calc->manip_->SetDOFValues(dof_vals);
  OR::Transform target = calc->pose_inv_.inverse(), cur = calc->link_->GetTransform();
  PlotAxes(env, cur, .05,  handles);
  PlotAxes(env, target, .05,  handles);
  handles.push_back(env.drawarrow(cur.trans, target.trans, .005, OR::Vector(1,0,1,1)));
}

VectorXd CartPoseConstraintCalculator::operator()(const VectorXd& dof_vals) const {
	manip_->SetDOFValues(toDblVec(dof_vals));
	OR::Vector curPosition = link_->GetTransform().trans;

	OR::Vector curError = curPosition - plane1_;

	Matrix<double, 1, 1> cost;
	cost[0] = curError.dot(normal_);

	return cost;
}


#if 0
struct CartPositionErrCalculator {
  Vector3d pt_world_;
  RobotAndDOFPtr manip_;
  OR::KinBody::LinkPtr link_;
  CartPositionErrCalculator(const Vector3d& pt_world, RobotAndDOFPtr manip, OR::KinBody::LinkPtr link) :
  pt_world_(pt_world),
  manip_(manip),
  link_(link)
  {}
  VectorXd operator()(const VectorXd& dof_vals) {
    manip_->SetDOFValues(toDblVec(dof_vals));
    OR::Transform newpose = link_->GetTransform();
    return pt_world_ - toVector3d(newpose.trans);
  }
};
#endif

MatrixXd CartVelJacCalculator::operator()(const VectorXd& dof_vals) const {
  int n_dof = manip_->GetDOF();
  MatrixXd out(6, 2*n_dof);
  manip_->SetDOFValues(toDblVec(dof_vals.topRows(n_dof)));
  OR::Transform pose0 = link_->GetTransform();
  MatrixXd jac0 = manip_->PositionJacobian(link_->GetIndex(), pose0.trans);
  manip_->SetDOFValues(toDblVec(dof_vals.bottomRows(n_dof)));
  OR::Transform pose1 = link_->GetTransform();
  MatrixXd jac1 = manip_->PositionJacobian(link_->GetIndex(), pose1.trans);
  out.block(0,0,3,n_dof) = -jac0;
  out.block(0,n_dof,3,n_dof) = jac1;
  out.block(3,0,3,n_dof) = jac0;
  out.block(3,n_dof,3,n_dof) = -jac1;
  return out;
}

VectorXd CartVelCalculator::operator()(const VectorXd& dof_vals) const {
  int n_dof = manip_->GetDOF();
  manip_->SetDOFValues(toDblVec(dof_vals.topRows(n_dof)));
  OR::Transform pose0 = link_->GetTransform();
  manip_->SetDOFValues(toDblVec(dof_vals.bottomRows(n_dof)));
  OR::Transform pose1 = link_->GetTransform();
  VectorXd out(6);
  out.topRows(3) = toVector3d(pose1.trans - pose0.trans - OR::Vector(limit_,limit_,limit_));
  out.bottomRows(3) = toVector3d( - pose1.trans + pose0.trans - OR::Vector(limit_, limit_, limit_));
  return out;
}


#if 0
CartVelConstraint::CartVelConstraint(const VarVector& step0vars, const VarVector& step1vars, RobotAndDOFPtr manip, KinBody::LinkPtr link, double distlimit) :
        ConstraintFromFunc(VectorOfVectorPtr(new CartVelCalculator(manip, link, distlimit)),
             MatrixOfVectorPtr(new CartVelJacCalculator(manip, link, distlimit)), concat(step0vars, step1vars), VectorXd::Ones(0), INEQ, "CartVel") 
{} // TODO coeffs
#endif

#if 0
struct UpErrorCalculator {
  Vector3d dir_local_;
  Vector3d goal_dir_world_;
  RobotAndDOFPtr manip_;
  OR::KinBody::LinkPtr link_;
  MatrixXd perp_basis_; // 2x3 matrix perpendicular to goal_dir_world
  UpErrorCalculator(const Vector3d& dir_local, const Vector3d& goal_dir_world, RobotAndDOFPtr manip, KinBody::LinkPtr link) :
    dir_local_(dir_local),
    goal_dir_world_(goal_dir_world),
    manip_(manip),
    link_(link)
  {
    Vector3d perp0 = goal_dir_world_.cross(Vector3d::Random()).normalized();
    Vector3d perp1 = goal_dir_world_.cross(perp0);
    perp_basis_.resize(2,3);
    perp_basis_.row(0) = perp0.transpose();
    perp_basis_.row(1) = perp1.transpose();
  }
  VectorXd operator()(const VectorXd& dof_vals) {
    manip_->SetDOFValues(toDblVec(dof_vals));
    OR::Transform newpose = link_->GetTransform();
    return perp_basis_*(toRot(newpose.rot) * dir_local_ - goal_dir_world_);
  }
};
#endif
}
