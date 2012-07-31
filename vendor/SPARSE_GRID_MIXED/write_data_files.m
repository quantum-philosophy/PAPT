function sparse_grid_mixed_write_tests ( )

%*****************************************************************************80
%
%% SPARSE_GRID_MIXED_WRITE_TESTS calls SPARSE_GRID_MIXED_WRITE with various arguments.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    05 March 2011
%
%  Author:
%
%    John Burkardt
%
%  Local Parameters:
%
%    Local, real TOL, a tolerance for point equality.
%    A value of sqrt ( eps ) is reasonable, and will allow the code to
%    consolidate points which are equal, or very nearly so.  A value of
%    -1.0, on the other hand, will force the code to use every point, regardless
%    of duplication.
%
  addpath ( '../sandia_rules' );

  tol = sqrt ( eps );

  fprintf ( 1, '\n' );
  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'SPARSE_GRID_MIXED_WRITE_TESTS\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '  Call SPARSE_GRID_MIXED_WRITE_TEST with various arguments.\n' );
  fprintf ( 1, '  All tests will use a point equality tolerance of %e\n', tol );

  for level_max = 0 : 10

    dim_num = 2;
%   level_max = 4;
    rule = [ 17, 17 ]';
    alpha = [ 0.0, 0.0 ]';
    beta = [ 0.0, 0.0 ]';
    file_name = sprintf ( 'ccn_d2_level%d', level_max );
    sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha, beta, tol, ...
      file_name );

  end
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'SPARSE_GRID_MIXED_WRITE_TESTS\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  rmpath ( '../sandia_rules' );

  return
end
