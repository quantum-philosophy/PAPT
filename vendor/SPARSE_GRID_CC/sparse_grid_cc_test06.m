function test06 ( dim_num, level_max )
 
%*****************************************************************************80
%
%% TEST06 creates a sparse Clenshaw-Curtis grid and writes it to a file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    11 August 2009
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer DIM_NUM, the spatial dimension.
%
%    Input, integer LEVEL_MAX, the level.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST06:\n' );
  fprintf ( 1, '  SPARSE_GRID_CC makes a sparse Clenshaw-Curtis grid.\n' );
  fprintf ( 1, '  Write the data to a set of quadrature files.\n' );
  
  fprintf ( 1, '\n' );
  fprintf ( 1, '  LEVEL_MAX = %d\n', level_max );
  fprintf ( 1, '  Spatial dimension DIM_NUM = %d\n', dim_num );
%
%  Determine the number of points.
%
  point_num = sparse_grid_cfn_size ( dim_num, level_max );
 
  r(1:dim_num,1) = -1.0;
  r(1:dim_num,2) = +1.0;
%
%  Compute the weights and points.
%
  [ w, x ] = sparse_grid_cc ( dim_num, level_max, point_num );
%
%  Write the data out.
%
  r_filename = sprintf ( 'cc_d%d_level%d_r.txt', dim_num, level_max );
  w_filename = sprintf ( 'cc_d%d_level%d_w.txt', dim_num, level_max );
  x_filename = sprintf ( 'cc_d%d_level%d_x.txt', dim_num, level_max );

  r8mat_write ( r_filename, dim_num, 2,         r );
  r8mat_write ( w_filename, 1,       point_num, w );
  r8mat_write ( x_filename, dim_num, point_num, x );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  R data written to "%s".\n', r_filename );
  fprintf ( 1, '  W data written to "%s".\n', w_filename );
  fprintf ( 1, '  X data written to "%s",\n', x_filename );

  return
end
