"""
#' Model predictions 
#' 
#' Function generates a data frame of model predictions for the typical value in the population,
#' individual predictions and data predictions.  The function can also be used to generate datasets
#' without predictions using the design specified in the arguments.
#' 
#' @param poped_db A PopED database created by \code{\link{create.poped.database}}.
#' @param models_to_use Which model numbers should we use? 
#' Model numbers are defined in \code{design} below using \code{model_switch}. For an explanation see \code{\link{create_design}}.
#' @param model_num_points How many extra observation rows should be created in the data frame for each group or individual 
#' per model.  If used then the points are placed evenly between \code{model_minxt} and \code{model_maxxt}. This option
#' is used by \code{\link{plot_model_prediction}} to simulate the response of the model on a finer grid then the defined design.
#' If \code{None} then only the input design is used.  Can be a single value or a vector the same length as the number of models.
#' @param model_minxt The minimum time value for extra observation rows indicated by \code{model_num_points}. 
#' A vector the same length as the number of models
#' @param model_maxxt The minimum time value for extra observation rows indicated by \code{model_num_points}. 
#' A vector the same length as the number of models
#' @param include_sample_times Should observations rows in the output data frame include the times indicated in the input design?
#' @param IPRED Should we simulate individual predictions?
#' @param DV should we simulate observations?
#' @param include_a Should we include the continuous design variables in the output?
#' @param include_x Should we include the discrete design variables in the output?
#' @param groups_to_use Which groups should we include in the output data frame?Allowed values are \code{"all"} or 
#' a vector of numbers indicating the groups to include, e.g. \code{c(1,3,6)}.
#' @param filename A filename that the data frame should be written to in comma separate value (csv) format.
#' @param dosing A list of lists that adds dosing records to the data frame (Each inner list corresponding to a group in the design). 
#' @param design A list that is passed as arguments to the function \code{\link{create_design}} to create a design object.  
#' @param model A list containing the model elements to use for the predictions
#' @param parameters A list of parameters to use in the model predictions.
#' @param predictions Should the resulting data frame have predictions?  Either \code{True} or \code{False} 
#' or \code{None} in which case the function decides based on other arguments.  
#' @param manipulation A list of one or more \code{\link[base]{expression}} arguments.  Each expression is 
#' evaluated using the code \code{for(i in 1:length(manipulation)){df = within(df,{eval(manipulation[[i]])})}}. 
#' Can be used to transform 
#' or create new columns in the resulting data frame. Note that these transformations are created after any model predictions occur,
#' so transformations in columns having to do with input to model predictions  will not affect the predictions.   
#' @param PI Compute prediction intervals for the data given the model.  Predictions are based on first-order approximations to 
#' the model variance and a normality assumption of that variance.
#' @param PI_conf_level The confidence level for the prediction interval computed.   
#' 
#' @return A dataframe containing a design and (potentially) simulated data with some dense grid of samples and/or 
#' based on the input design.
#'  
#' @family evaluate_design
#' @family Simulation
#' 
#' @example tests/testthat/examples_fcn_doc/examples_model_prediction.R
#' 
#' @export
# @importFrom mvtnorm rmvnorm
# @importFrom dplyr rbind_list
Author: Caiya Zhang, Yuchen Zheng
"""


from re import M
from project.param_choose import param_choose
from project.poped_choose import poped_choose
from project.create_design import create_design


def model_prediction(poped_db=None,
                     design={},
                     model={},
                     parameters={},
                     IPRED=False,
                     DV=False,
                     dosing=None,
                     predictions=None,
                     filename=None,
                     models_to_use="all",
                     model_num_points=None,
                     model_minxt=None,
                     model_maxxt=None,
                     include_sample_times=True,
                     groups_to_use="all",
                     include_a=True,
                     include_x=True,
                     manipulation=None,
                     PI=False,
                     PI_conf_level=0.95):
    if len(design) == 0:
        design = {"xt": poped_db["design"]["xt"],
                  "groupsize": poped_db["design"]["groupsize"],
                  "m": poped_db["design"]["m"],
                  "x": poped_db["design"]["x"],
                  "a": poped_db["design"]["a"],
                  "ni": poped_db["design"]["ni"],
                  "model_switch": poped_db["design"]["model_switch"]}
    if len(model) == 0:
        model = {"fg_pointer": poped_db["model"]["fg_pointer"],
                 "ff_pointer": poped_db["model"]["ff_pointer"],
                 "ferror_pointer": poped_db["model"]["ferror_pointer"]}
    if len(parameters) == 0:
        parameters = {"docc": poped_db["parameters"]["docc"],
                      "d": poped_db["parameters"]["d"],
                      "bpop": poped_db["parameters"]["bpop"],
                      "covd": poped_db["parameters"]["covd"],
                      "covdocc": poped_db["parameters"]["covdocc"],
                      "sigma": poped_db["parameters"]["sigma"]}
    if predictions is None:
        predictions = False
        if poped_db is not None or len(parameters) != 0 and len(model) != 0 and len(design) != 0:
            predictions = True
    
    if poped_db is None and len(design) == 0:
        raise Exception("Either 'poped_db' or 'design' need to be defined")
    
    design = create_design(design)

    NumOcc = poped_db["parameters"]["NumOcc"]

    if DV:
        IPRED = True

    # with design
    maxxt = param_choose(poped_db, design["xt"], 0, argv=["design_space", "maxxt"]) # Matrix defining the max value for each sample
    minxt = param_choose(poped_db, design["xt"], 0, argv=["design_space", "minxt"]) # Matrix defining the min value for each sample

    # size checking
    if dosing is not None:
        if len(dosing) != design["m"]:
            if len(dosing) == 1:
                dosing = dosing * design["m"]
            else:
                raise Exception("dosing argument does not have the right dimensions. Must be 1 list or a list of lists the size of the number of groups")
    
    if predictions:
        docc_size = 0