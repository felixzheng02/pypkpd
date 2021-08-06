"""
Plot model predictions 

Function plots model predictions for the typical value in the population,
individual predictions and data predictions_

@inheritParams RS_opt
@inheritParams model_prediction
@param separate_groups Should there be separate plots for each group_
@param sample_times Should sample times be shown on the plots_
@param sample_times_IPRED Should sample times be shown based on the IPRED y-values_
@param sample_times_DV Should sample times be shown based on the DV y-values_
@param PRED Should a PRED line be drawn_
@param IPRED_lines Should IPRED lines be drawn?
@param alpha_IPRED_lines What should the transparency for the IPRED_lines be?
@param alpha_IPRED What should the transparency of the IPRED CI?
@param sample_times_size What should the size of the sample_times be?
@param alpha_DV What should the transparency of the DV CI?
@param DV_lines Should DV lines be drawn?
@param DV_points Should DV points be drawn?
@param alpha_DV_lines What should the transparency for the DV_lines be?
@param alpha_DV_points What should the transparency for the DV_points be?
@param sample_times_DV_points TRUE or FALSE_
@param sample_times_DV_lines TRUE or FALSE_
@param alpha_sample_times_DV_points What should the transparency for the sample_times_DV_points be?
@param alpha_sample_times_DV_lines What should the transparency for the sample_times_DV_lines be?
@param y_lab The label of the y-axis_
@param facet_scales Can be "free", "fixed", "free_x" or "free_y"
@param facet_label_names TRUE or FALSE
@param IPRED_lines_pctls Should lines be drawn at the chosen percentiles of the IPRED values?  
@param groupsize_sim How many individuals per group  should be 
  simulated when DV=TRUE or IPRED=TRUE to create prediction intervals?
@param model_names A vector of names of the response model/s (the length of the 
vector should be equal to the number of response models)_ It is Null by default_
@param PI Plot prediction intervals for the expected data given the model_  
Predictions are based on first-order approximations to 
the model variance and a normality assumption of that variance_  As such these computations are 
more approximate than using \code{DV=T} and \code{groupsize_sim = some large number}_
# @param PI_fill The color of the PI_
@param PI_alpha The transparency of the PI_
@param DV_mean_sd Plot the mean and standard deviation of simulated observations_ 
@param ___ Additional arguments passed to the \code{\link{model_prediction}} function_

@return A \link[ggplot2]{ggplot} object_  If you would like to further edit this plot don't 
forget to load the ggplot2 library using \code{library(ggplot2)}_

@family evaluate_design
@family Simulation
@family Graphics

@seealso \code{\link{model_prediction}}

@example tests/testthat/examples_fcn_doc/examples_plot_model_prediction_R

@export
@import ggplot2
# @import Hmisc
"""


from project.model_prediction import model_prediction


def plot_model_prediction(poped_db,
						model_num_points=100,
						groupsize_sim = 100,
						separate_groups=False, 
						sample_times=True, 
						sample_times_IPRED=False,
						sample_times_DV=False,
						PRED=True,
						IPRED=False,
						IPRED_lines=False,
						IPRED_lines_pctls=False,
						alpha_IPRED_lines=0.1,
						alpha_IPRED=0.3,
						sample_times_size=4,
						DV=False,
						alpha_DV=0.3,
						DV_lines=False,
						DV_points=False,
						alpha_DV_lines=0.3,
						alpha_DV_points=0.3,
						sample_times_DV_points=False,
						sample_times_DV_lines=False,
						alpha_sample_times_DV_points=0.3,
						alpha_sample_times_DV_lines=0.3,
						y_lab="Model Predictions",
						facet_scales="fixed", # could be "free", "fixed", "free_x" or "free_y"
						facet_label_names = True, 
						model_names=None,
						DV_mean_sd=False,
						PI=False,
						PI_alpha=0.3,
						*argv):
	PI_u = None
	PI_l = None
	df = model_prediction(poped_db,
						  argv,
						  model_num_points=model_num_points,
						  PI=PI
						  )
	if model_names is not None:
		...
		# levels(df.2$Model) <- model.names
	
	if IPRED or IPRED_lines or DV or IPRED_lines_pctls or sample_times_DV or sample_times_DV_points or sample_times_DV_lines or DV_mean_sd:
		dv_val = False
		if DV or sample_times_DV or sample_times_DV_points or sample_times_DV_lines or DV_mean_sd:
			dv_val = True
		poped_db_tmp = poped_db
		poped_db_tmp["design"]["groupsize"] = poped_db["design"]["groupsize"] * 0 + groupsize_sim
		df_ipred = model_prediction(poped_db_tmp,
									argv,
									model_num_points=model_num_points,
									IPRED=True,
									DV=dv_val)
		if model_names is not None:
			...
			# levels(df.ipred$Model) <- model.names
		if sample_times_IPRED or sample_times_DV or sample_times_DV_points or sample_times_DV_lines:
			# df.ipred.samples <- df.ipred[df.ipred$Time %in% poped.db$design$xt,]
			if model_names is not None:
				...
				# levels(df.ipred.samples$Model) <- model.names