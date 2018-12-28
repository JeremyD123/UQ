function [ grid_weight, grid_point ] = sparse_grid_laguerre ( dim_num, ...
  level_max, point_num )
  
%*****************************************************************************80
%
%% SPARSE_GRID_LAGUERRE computes a sparse grid of Gauss-Laguerre points. 
% 
%  Discussion: 
% 
%    The quadrature rule is associated with a sparse grid derived from 
%    a Smolyak construction using a 1D Gauss-Laguerre quadrature rule.  
% 
%    The user specifies: 
%    * the spatial dimension of the quadrature region, 
%    * the level that defines the Smolyak grid. 
% 
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    05 July 2008
% 
%  Author: 
% 
%    John Burkardt 
%
%  Reference:
%
%    Fabio Nobile, Raul Tempone, Clayton Webster,
%    A Sparse Grid Stochastic Collocation Method for Partial Differential
%    Equations with Random Input Data,
%    SIAM Journal on Numerical Analysis,
%    Volume 46, Number 5, 2008, pages 2309-2345.
%
%  Parameters: 
% 
%    Input, integer DIM_NUM, the spatial dimension. 
% 
%    Input, integer LEVEL_MAX, controls the size of the 
%    sparse grid. 
% 
%    Input, integer POINT_NUM, the number of points in the grid, 
%    as determined by SPARSE_GRID_LAGUERRE_SIZE. 
% 
%    Output, real GRID_WEIGHT(POINT_NUM), the weights. 
% 
%    Output, real GRID_POINT(DIM_NUM,POINT_NUM), the points. 
% 
  grid_weight(1:point_num) = 0.0;
%
%  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
%
  point_num2 = 0;

  level_min = max ( 0, level_max + 1 - dim_num );
  
  for level = level_min : level_max
%
%  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
%  that adds up to LEVEL.
%
    level_1d = [];
    more = 0;
    h = 0;
    t = 0;

    while ( 1 )

      [ level_1d, more, h, t ] = comp_next ( level, dim_num, level_1d, more, h, t );
%
%  Transform each 1D level to a corresponding 1D order.
%  The relationship is the same as for other OPEN rules.
%
      order_1d = level_to_order_open ( dim_num, level_1d );

      grid_base2(1:dim_num) = order_1d(1:dim_num);
%
%  The product of the 1D orders gives us the number of points in this subgrid.
%
      order_nd = prod ( order_1d(1:dim_num) );
%
%  Compute the weights for this product grid.
%
      grid_weight2 = product_weight_laguerre ( dim_num, order_1d, order_nd );
%
%  Now determine the coefficient of the weight.
%
      coeff = (-1)^( level_max - level ) ...
        * i4_choose ( dim_num - 1, level_max - level );
%
%  The inner (hidden) loop generates all points corresponding to given grid.
%  The grid indices will be between 1 and ORDER_1D(DIM).
%
      grid_index2 = multigrid_index_one ( dim_num, order_1d, order_nd );

      for point = 1 : order_nd

        point_num2 = point_num2 + 1;

        grid_point(1:dim_num,point_num2) = laguerre_abscissa ( dim_num, 1, ...
          grid_index2(1:dim_num,point), grid_base2(1:dim_num) );

        grid_weight(point_num2) = coeff * grid_weight2(point);

      end

      if ( ~more )
        break
      end

    end

  end

  return
end
