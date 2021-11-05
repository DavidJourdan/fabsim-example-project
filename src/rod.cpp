#include <fsim/ElasticRod.h>
#include <fsim/RodCollection.h>
#include <fsim/util/typedefs.h>
#include <optim/NewtonSolver.h>
#include <polyscope/curve_network.h>

int main(int argc, char *argv[]) 
{
  using namespace Eigen;

  // initialize geometry
  fsim::Mat3<double> V(94, 3);
  for (int i = 0; i < 94; ++i)
    // V.row(i) << sin(i / 15.), -cos(i / 15.), 0;
    V.row(i) << i / 94., 0, 0; // uncomment this line for a cantilever test

  // parameters of the rod model, cross-section can be circular or rectangular
  const double thickness = 1;
  const double width = 7;
  const double young_modulus = 10;
  const double mass = 1;
  fsim::RodParams params = {thickness, width, young_modulus, mass, fsim::CrossSection::Square};
  
  // declare ElasticRod object
  fsim::ElasticRod rod(V, Vector3d::UnitZ(), params);

  // declare NewtonSolver object
  optim::NewtonSolver<double> solver;
  // specify fixed degrees of freedom (here: first two vertices and first frame rotation are fixed)
  solver.options.fixed_dofs = {0 * 3 + 0, 0 * 3 + 1, 0 * 3 + 2, 1 * 3 + 0,
                               1 * 3 + 1, 1 * 3 + 2, 3 * 94};
  solver.options.threshold = 1e-6; // specify how small the gradient's norm has to be
  solver.options.update_fct = [&](const Ref<const VectorXd> X) {
    rod.updateProperties(X); // updates twist angles as per [Bergou et al. 2010] (section 6)
  };

  // display the curve and its material frames
  polyscope::registerCurveNetworkLine("rod", V);
  fsim::Mat3<double> D1, D2;
  rod.getReferenceDirectors(D1, D2);
  polyscope::getCurveNetwork("rod")
      ->addEdgeVectorQuantity("D1", D1)
      ->setEnabled(true);
  polyscope::getCurveNetwork("rod")
      ->addEdgeVectorQuantity("D2", D2)
      ->setEnabled(true);
  polyscope::view::upDir = polyscope::view::UpDir::ZUp;
  polyscope::options::groundPlaneHeightFactor = 0.4;
  polyscope::init();

  polyscope::state::userCallback = [&]() 
  {
    ImGui::PushItemWidth(100);
    if(ImGui::InputDouble("Thickness", &params.thickness, 0, 0, "%.1f"))
      rod.setParams(params);
    
    if(ImGui::InputDouble("Width", &params.width, 0, 0, "%.1f"))
      rod.setParams(params);
    
    if(ImGui::InputDouble("Young's modulus", &params.E, 0, 0, "%.1f"))
      rod.setParams(params);
    
    if(ImGui::InputDouble("Mass", &params.mass, 0, 0, "%.1f"))
      rod.setParams(params);

    if(ImGui::Button("Solve")) 
    {
      // prepare input variables for the Newton solver
      VectorXd var = VectorXd::Zero(V.size() + V.rows() - 1);
      var.head(V.size()) = Map<VectorXd>(V.data(), V.size());
      rod.updateProperties(var); // necessary to do so if it's not the first solve
      
      // Newton's method: finds a local minimum of the energy (Fval = energy value, Optimality = gradient's norm)
      solver.solve(rod, var);

      // Display the result of the optimization
      polyscope::getCurveNetwork("rod")->updateNodePositions(
          Map<fsim::Mat3<double>>(solver.var().data(), V.rows(), 3));

      rod.getRotatedDirectors(solver.var().tail(V.rows() - 1), D1, D2);
      polyscope::getCurveNetwork("rod")->addEdgeVectorQuantity("D1", D1);
      polyscope::getCurveNetwork("rod")->addEdgeVectorQuantity("D2", D2);
    }
  };
  polyscope::show();
}
