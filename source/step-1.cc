/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * based on deal.II step-1
 */


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <cmath>
#include <fstream>
#include <iostream>

using namespace dealii;


//! Generate a hypercube, and output it as an svg file.
void first_grid(Triangulation<2> &triangulation)
{
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(4);

  std::ofstream out("grid-1.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to grid-1.svg" << std::endl;
}


//! Generate a locally refined hyper_shell, and output it as an svg file.
void second_grid(Triangulation<2> &triangulation)
{
  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 10);


  // triangulation.reset_manifold(0);

  for (unsigned int step = 0; step < 5; ++step)
    {
      for (auto &cell : triangulation.active_cell_iterators())
        {
          for (const auto v : cell->vertex_indices())
            {
              const double distance_from_center =
                center.distance(cell->vertex(v));

              if (std::fabs(distance_from_center - inner_radius) <=
                  1e-6 * inner_radius)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }

      triangulation.execute_coarsening_and_refinement();
    }


  std::ofstream out("grid-2.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);

  std::cout << "Grid written to grid-2.svg" << std::endl;
}


//! Create an L-shaped domain with one global refinement, and write it on
// `third_grid.vtk`.  Refine the L-shaped mesh adaptively around the re-entrant
// corner three times (after the global refinement you already did), but with a
// twist: refine all cells with the distance between the center of the cell and
// re-entrant corner is smaller than 1/3.
void third_grid(Triangulation<2> &tria)
{
  // Insert code here
  // last year code in comments
  /* Triangulation<2> tr1,tr2,tr3,tr_final;
   const Point<2> p1(0, 0);
   const Point<2> p2(1,2);
   const Point<2> p3(2,1);
   const Point<2> p4(1,0);
   const Point<2> p5(1,1);
   const Point<2> p6(0,1);
   GridGenerator::hyper_rectangle(tr1,p1,p5);
   GridGenerator::hyper_rectangle(tr2,p6,p2);
   GridGenerator::hyper_rectangle(tr3,p4,p3);
   //they cannot be refined!
   GridGenerator::merge_triangulations ({&tr1,&tr2,&tr3},tr_final);*/
  GridGenerator::hyper_L(tria);
  tria.refine_global(1);
  const Point<2> p5(0, 0);
  for (unsigned int step = 0; step < 3; ++step)
    {
      // Active cells are those that are not further refined
      // we need to mark cells for refinement
      for (auto &cell : tria.active_cell_iterators())
        {
          const double distance_from_corner = cell->center().distance({0, 0});
          // choose whatever refinement condition
          if (distance_from_corner < 2.0 / 3.0)
            {
              cell->set_refine_flag();
            }
        }
      // refine global calls this function too
      tria.execute_coarsening_and_refinement();
    } // for steps loop

  // tr_final.refine_global(2); //we want a nice picture :)

  std::ofstream file_var("grid-3.vtk");
  GridOut       grid_out;
  grid_out.write_vtk(tria, file_var);
  std::cout << "Grid written to grid-3.vtk" << std::endl;
}

//! Returns a tuple with number of levels, number of cells, number of active
// cells. Test this with all of  your meshes.
std::tuple<unsigned int, unsigned int, unsigned int>
get_info(const Triangulation<2> &tria)
{
  // Insert code here
  return std::make_tuple(tria.n_levels(),
                         tria.n_cells(),
                         tria.n_active_cells());
}

void
torus_grid()
{
  Triangulation<2, 3> tria;
  GridGenerator::torus(tria, 2, 1);
  tria.refine_global(2);

  std::ofstream out("grid-torus.vtk");
  GridOut       grid_out;
  grid_out.write_vtk(tria, out);
  std::cout << "Grid written to grid-torus.vtk" << std::endl;
}


int
main()
{
  Triangulation<2> triangulation;
  first_grid(triangulation);
  triangulation.clear();
  second_grid(triangulation);
  triangulation.clear();
  third_grid(triangulation);

  torus_grid();
}
