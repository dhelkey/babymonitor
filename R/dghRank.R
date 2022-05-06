dghRank = function (y_vec, mean_vec, var_vec, ind_mat)
{
	#' Calculate Standardized Values
	#'
    #' Optimized by Lucy Greenberg
	#' @param y_vec <nx1> Individual observation vector
	#' @param mean_vec <nx1> Predicted individual values
	#' @param var_vec <nx1> Prediction variance
	#' @param ind_mat <nxp> Institution membership matrix
    if (is.null(ind_mat)) return(NULL)
      n_inv = diag(1/Matrix::colSums(ind_mat))
      n_inv = Matrix(n_inv,sparse = TRUE)
          p = dim(ind_mat)[2]
          t_ind_mat<-Matrix::t(ind_mat)
          e_vec = matrix(n_inv %*% (t_ind_mat %*% mean_vec))
          o_vec = matrix(n_inv %*% (t_ind_mat %*% y_vec))
          d_vec = o_vec - e_vec
          s_vec = matrix(n_inv %*% sqrt(t_ind_mat %*% var_vec))
          return(list(E = e_vec,
                      O = o_vec,
                      D = d_vec,
                      S = s_vec,
                      Z = d_vec/s_vec
                      ))

}
