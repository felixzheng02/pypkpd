"""
## Compute the Bayesian Fisher information Matrix 
## 
## Computation of the Bayesian Fisher information Matrix for 
## individual parameters of a population model based on 
## Maximum A Posteriori (MAP) estimation of the empirical Bayes estimates (EBEs) in a population model 
##
## @param poped_db A PopED database
## @param num_sim_ids If \code{use_mc=TRUE}, how many individuals should be
##   simulated to make the computations.
## @param use_mc Should the calculation be based on monte-carlo simulations. If
##   not then then a first order approximation is used
## @param use_purrr If \code{use_mc=TRUE} then should the method use the package
##   purrr in calculations?  This may speed up computations (potentially).
## @param shrink_mat Should the shrinkage Matrix be returned.  Calculated as the
## inverse of the  Bayesian Fisher information Matrix times the inverse of the 
## omega Matrix (variance Matrix of the between-subject variability).
## @return The Bayesian Fisher information Matrix for each design group 
## @export
##
## @references \enumerate{ 
##   \item Combes, F. P., Retout, S.,
##   Frey, N., & Mentre, F. (2013). Prediction of shrinkage of individual
##   parameters using the Bayesian information Matrix in non-linear mixed effect
##   models with evaluation in pharmacokinetics. Pharmaceutical Research, 30(9),
##   2355-67. \doi{10.1007/s11095-013-1079-3}. 
##   \item Hennig, S., Nyberg, J., Fanta, S., Backman, J.
##   T., Hoppu, K., Hooker, A. C., & Karlsson, M. O. (2012). Application of the
##   optimal design approach to improve a pretransplant drug dose finding design
##   for ciclosporin. Journal of Clinical Pharmacology, 52(3), 347-360.
##   \doi{10.1177/0091270010397731}. 
##   }

## Author: Caiya Zhang, Yuchen Zheng
"""


from project.create_poped_database import create_poped_database
from project.find_largest_index import find_largest_index
from project.capture_output import capture_output

def evaluate_fim_map(poped_db,
                     use_mc=False,
                     num_sim_ids=1000,
                     use_purrr=False,
                     shrink_mat=False):
    group = None

    # tranform random effects to fixed effects
    tmp_fg = poped_db["model"]["fg_pointer"]
    if type(tmp_fg) is str:
        tmp_fg = eval(tmp_fg)
    largest_bpop = find_largest_index(tmp_fg, lab="b")
    largest_eps = find_largest_index(poped_db["model"]["ferror_pointer"], lab="epsi", mat=True, mat_row=False)
    txt = capture_output()


def replace_fun(txt):
    ...