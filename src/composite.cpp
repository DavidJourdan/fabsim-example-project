#include <fsim/CompositeModel.h>
#include <fsim/ElasticRod.h>
#include <fsim/ElasticMembrane.h>
#include <fsim/util/io.h>
#include <fsim/util/typedefs.h>
#include <optim/NewtonSolver.h>
#include <polyscope/curve_network.h>
#include <polyscope/surface_mesh.h>
#include <algorithm>

int main(int argc, char *argv[]) {
  using namespace Eigen;

  // load geometry from OFF mesh file
  fsim::Mat3<double> V;
  fsim::Mat3<int> F;
  fsim::readOFF("../data/mesh.off", V, F);

  // parameters of the membrane model
  const double young_modulus = 10;
  const double thickness = 0.5;
  const double poisson_ratio = 0.3;
  double stretch_factor = 1.2;
  double mass = 0;

  RowVector3d barycenter = V.colwise().sum() / V.rows();
  VectorXi indices;
  for(int k = 0; k < 4; ++k)
  {
    // find vertices in edge k of the mesh
    std::vector<int> indicesInEdge;
    for(int i = 0; i < V.rows(); ++i)
      if(fsim::point_in_segment(V.row(i), V.row(k), V.row((k + 1) % 4)))
      {
        indicesInEdge.push_back(i);
      }
    // sort vertices
    std::sort(indicesInEdge.begin(), indicesInEdge.end(), [&](int a, int b) { 
      Vector3d orientation = (V.row(a) - barycenter).cross(V.row(b) - barycenter);
      return orientation(2) > 0; 
    });

    int n = indices.size();
    indices.conservativeResize(n + indicesInEdge.size() - 1);
    indices.segment(n, indicesInEdge.size() - 1) = Map<VectorXi>(indicesInEdge.data(), indicesInEdge.size() - 1);
  }

  fsim::RodParams params = {1, 1.5, 1000 * young_modulus, mass, fsim::CrossSection::Circle, true};

  fsim::CompositeModel composite(
    fsim::StVKMembrane(V / stretch_factor, F, thickness, young_modulus, poisson_ratio, mass), 
    fsim::ElasticRod(V, indices, Vector3d::UnitZ(), params));

  // declare NewtonSolver object
  optim::NewtonSolver<double> solver;
  // // specify fixed degrees of freedom (here the 4 corners of the mesh are fixed)
  // solver.options.fixed_dofs = {0 * 3 + 0, 0 * 3 + 1, 0 * 3 + 2, 1 * 3 + 0,
  //                              1 * 3 + 1, 1 * 3 + 2, 2 * 3 + 0, 2 * 3 + 1,
  //                              2 * 3 + 2, 3 * 3 + 0, 3 * 3 + 1, 3 * 3 + 2};
  solver.options.threshold = 1e-6; // specify how small the gradient's norm has to be
  solver.options.update_fct = [&](const Ref<const VectorXd> X) {
    // updates twist angles as per [Bergou et al. 2010] (section 6)
    composite.getModel<1>().updateProperties(X); 
  };

  // display the curve network
  fsim::Mat3<double> R(indices.size() - 1, 3);
  fsim::Mat2<int> E(indices.size() - 1, 2);
  for(int j = 0; j < indices.size() - 1; ++j)
  {
    R.row(j) = V.row(indices(j));
    E.row(j) << j, j + 1;
  }
  E(E.rows() - 1, 1) = 0;
  polyscope::registerCurveNetwork("rod", R, E);
  
  // display the mesh
  polyscope::registerSurfaceMesh("mesh", V, F)
      ->setEdgeWidth(1)
      ->setEdgeColor({0.1, 0.1, 0.1});
  polyscope::view::upDir = polyscope::view::UpDir::ZUp;
  polyscope::options::groundPlaneHeightFactor = 0.4;
  polyscope::init();

  polyscope::state::userCallback = [&]() 
  {
    // ImGui::PushItemWidth(100);
    // if(ImGui::InputDouble("Stretch factor", &stretch_factor, 0, 0, "%.1f"))
    //   membrane = fsim::StVKMembrane(V.leftCols(2) / stretch_factor, F, thickness, young_modulus, poisson_ratio, mass);
    
    // if(ImGui::InputDouble("Mass", &mass, 0, 0, "%.1f"))
    //   membrane.setMass(mass);

    if(ImGui::Button("Solve")) 
    {
      // prepare input variables for the Newton solver
      VectorXd var = VectorXd::Zero(V.size() + V.rows());
      var.head(V.size()) = Map<VectorXd>(V.data(), V.size());

      for(int i = 0; i < V.rows(); ++i)
      {
        // add noise in the Z direction to force the rod out of plane
        std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<double> dis(-0.1, 0.1);
        var(3 * i + 2) = dis(gen);
      }

      composite.getModel<1>().updateProperties(var); // necessary to do so if it's not the first solve

      // Newton's method: finds a local minimum of the energy (Fval = energy value, Optimality = gradient's norm)
      var = solver.solve(composite, var);

      // Display the result of the optimization
      polyscope::getSurfaceMesh("mesh")->updateVertexPositions(
          Map<fsim::Mat3<double>>(var.data(), V.rows(), 3));

      for(int j = 0; j < indices.size() - 1; ++j)
        R.row(j) = var.segment<3>(3 * indices(j));
      polyscope::getCurveNetwork("rod")->updateNodePositions(R);
    }
  };
  polyscope::show();
}
