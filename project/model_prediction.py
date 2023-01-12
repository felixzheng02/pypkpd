"""
## Model predictions 
## 
## Function generates a data frame of model predictions for the typical value in the population,
## individual predictions and data predictions.  The function can also be used to generate datasets
## without predictions using the design specified in the arguments.
## 
## @param pypkpd_db A PopED database created by \code{\link{create.poped.database}}.
## @param models_to_use Which model numbers should we use? 
## Model numbers are defined in \code{design} below using \code{model_switch}. For an explanation see \code{\link{create_design}}.
## @param model_num_points How many extra observation rows should be created in the data frame for each group or individual 
## per model.  If used then the points are placed evenly between \code{model_minxt} and \code{model_maxxt}. This option
## is used by \code{\link{plot_model_prediction}} to simulate the response of the model on a finer grid then the defined design.
## If \code{None} then only the input design is used.  Can be a single value or a vector the same length as the number of models.
## @param model_minxt The minimum time value for extra observation rows indicated by \code{model_num_points}. 
## A vector the same length as the number of models
## @param model_maxxt The minimum time value for extra observation rows indicated by \code{model_num_points}. 
## A vector the same length as the number of models
## @param include_sample_times Should observations rows in the output data frame include the times indicated in the input design?
## @param IPRED Should we simulate individual predictions?
## @param DV should we simulate observations?
## @param include_a Should we include the continuous design variables in the output?
## @param include_x Should we include the discrete design variables in the output?
## @param groups_to_use Which groups should we include in the output data frame?Allowed values are \code{"all"} or 
## a vector of numbers indicating the groups to include, e.g. \code{c(1,3,6)}.
## @param filename A filename that the data frame should be written to in comma separate value (csv) format.
## @param dosing A list of lists that adds dosing records to the data frame (Each inner list corresponding to a group in the design). 
## @param design A list that is passed as arguments to the function \code{\link{create_design}} to create a design object.  
## @param model A list containing the model elements to use for the predictions
## @param parameters A list of parameters to use in the model predictions.
## @param predictions Should the resulting data frame have predictions?  Either \code{True} or \code{False} 
## or \code{None} in which case the function decides based on other arguments.  
## @param manipulation A list of one or more \code{\link[base]{expression}} arguments.  Each expression is 
## evaluated using the code \code{for(i in 1:length(manipulation)){df = within(df,{eval(manipulation[[i]])})}}. 
## Can be used to transform 
## or create new columns in the resulting data frame. Note that these transformations are created after any model predictions occur,
## so transformations in columns having to do with input to model predictions  will not affect the predictions.   
## @param PI Compute prediction intervals for the data given the model.  Predictions are based on first-order approximations to 
## the model variance and a normality assumption of that variance.
## @param PI_conf_level The confidence level for the prediction interval computed.   
## 
## @return A dataframe containing a design and (potentially) simulated data with some dense grid of samples and/or 
## based on the input design.
##  
## @family evaluate_design
## @family Simulation
## 
## @example tests/testthat/examples_fcn_doc/examples_model_prediction.R
## 
## @export
# @importFrom mvtnorm rmvnorm
# @importFrom dplyr rbind_list
Author: Caiya Zhang, Yuchen Zheng
"""


from matpy.matrix import Matrix
import numpy as np
import pandas as pd
from re import M
from project.param_choose import param_choose
from project.poped_choose import poped_choose
from.pypkpd_choose import pypkpd_choose
from project.create_design import create_design
from project.util import get_dict_value
from project.length import length
from project.data import data
from project.size import size
from project.feval import feval
from project.getfulld import getfulld


def model_prediction(pypkpd_db = None,
                     design = None,
                     model = None,
                     parameters = None,
                     IPRED = False,
                     DV = False,
                     dosing = None,
                     predictions = None,
                     filename = None,
                     models_to_use = "all",
                     model_num_points = None,
                     model_minxt = None,
                     model_maxxt = None,
                     include_sample_times = True,
                     groups_to_use = "all",
                     include_a = True,
                     include_x = True,
                     manipulation = None,
                     PI = False,
                     PI_conf_level = 0.95):

    design = param_choose(design, 
                            {"xt": get_dict_value(pypkpd_db, "design", "xt"),
                            "groupsize": get_dict_value(pypkpd_db, "design", "groupsize"),
                            "m": get_dict_value(pypkpd_db, "design", "m"),
                            "x": get_dict_value(pypkpd_db, "design", "x"),
                            "a": get_dict_value(pypkpd_db, "design", "a"),
                            "ni": get_dict_value(pypkpd_db, "design", "ni"),
                            "model_switch": get_dict_value(pypkpd_db, "design", "model_switch")}, 
                            Exception)
    model = param_choose(model, 
                            {"fg_pointer": get_dict_value(pypkpd_db, "model", "fg_pointer"),
                            "ff_pointer": get_dict_value(pypkpd_db, "model", "fg_pointer"),
                            "ferror_pointer": get_dict_value(pypkpd_db, "model", "ferror_pointer")}, 
                            Exception)
    parameters = param_choose(parameters,
                                {"docc": get_dict_value(pypkpd_db, "parameters", "docc"),
                                "d": get_dict_value(pypkpd_db, "parameters", "d"),
                                "bpop": get_dict_value(pypkpd_db, "parameters", "bpop"),
                                "covd": get_dict_value(pypkpd_db, "parameters", "covd"),
                                "covdocc": get_dict_value(pypkpd_db, "parameters", "covdocc"),
                                "sigma": get_dict_value(pypkpd_db, "parameters", "sigma")},
                                Exception)
    
    if predictions is None:
        predictions = False
        if pypkpd_db is not None or len(parameters) != 0 and len(model) != 0 and len(design) != 0:
            predictions = True
    
    if pypkpd_db is None and len(design) == 0:
        raise Exception("Either 'pypkpd_db' or 'design' need to be defined")
    
    design = create_design(xt = get_dict_value(design, "xt"),
                            groupsize = get_dict_value(design, "groupsize"),
                            m = get_dict_value(design, "m"),
                            x = get_dict_value(design, "x"),
                            a = get_dict_value(design, "a"),
                            ni = get_dict_value(design, "ni"),
                            model_switch = get_dict_value(design, "model_switch"))

    NumOcc = get_dict_value(pypkpd_db, "parameters", "NumOcc")

    if DV:
        IPRED = True

    # with design
    xt = get_dict_value(design, "xt")
    m = data(get_dict_value(design, "m"))
    ni = get_dict_value(design, "ni")
    model_switch = get_dict_value(design, "model_switch")
    a = get_dict_value(design, "a")
    x = get_dict_value(design, "x")
    groupsize = get_dict_value(design, "groupsize")

    maxxt = pypkpd_choose(get_dict_value(pypkpd_db, "design_space", "maxxt"), xt) # Matrix defining the max value for each sample
    minxt = pypkpd_choose(get_dict_value(pypkpd_db, "design_space", "minxt"), xt) # Matrix defining the min value for each sample

    # size checking
    if dosing is not None:
        if length(dosing) != m:
            if length(dosing) == 1:
                dosing = dosing * m
            else:
                raise Exception("dosing argument does not have the right dimensions. Must be 1 list or a list of lists the size of the number of groups")
    
    if predictions:
        docc_size = 0
        # omitted
        d_size = 0
        if length(data(get_dict_value(parameters, "d"))[:, 1]) != 0:
            d_size = size(data(get_dict_value(parameters, "d"))[:, 1])[0]
    
    used_times = np.zeros(size(xt))
    for i in range(0, size(xt)):
        used_times[i, 0:data(ni)[0, i]] = 1
    if groups_to_use == "all":
        groups_to_use = list(range(0, size(xt)[0]))
    if models_to_use == "all":
        models_to_use = np.unique(data(model_switch)[used_times == 1])
    
    id_num_start = 1
    for i in range(0, length(groups_to_use)):
        if a is None:
            a_i = np.zeros([0, 1])
        else:
            a_i = data(a)[data(groups_to_use)[i], :]
        if x is None:
            x_i = np.zeros([0, 1])
        else:
            x_i = data(x)[data(groups_to_use)[i], :]
        num_ids = groupsize[data(groups_to_use)[i]]
        if model_num_points is None:
            xt_i = data(xt)[data(groups_to_use)[i], 0:data(ni)[data(groups_to_use)[i]]]
            model_switch_i = data(model_switch)[data(groups_to_use)[i], 0:data(ni)[data(groups_to_use)[i]]]
            # omitted
        else:
            if length(models_to_use) > 1 and length(model_num_points) == 1:
                model_num_points = list(model_num_points) * length(models_to_use)
            for j in list(models_to_use):

                if model_minxt is None:
                    minv = np.min(data(minxt)[data(model_switch) == j])
                else:
                    if length(models_to_use) > 1 and length(model_minxt) == 1:
                        model_minxt = list(model_minxt) * length(models_to_use)
                    minv = list(model_minxt)[j]

                if model_maxxt is None:
                    maxv = np.max(data(maxxt)[data(model_switch) == j])
                else:
                    if length(models_to_use) > 1 and length(model_maxxt) == 1:
                        model_maxxt = list(model_maxxt) * length(models_to_use)
                    maxv = list(model_maxxt)[j]
                tmp_num_pts = model_num_points[j]
                if length(model_num_points) < j:
                    tmp_num_pts = model_num_points[0]
                xt_i = np.concatenate(data(xt_i), np.array(range(minv, maxv)), axis=0)
                model_switch_i = np.concatenate(data(model_switch_i), j*np.ones(1, tmp_num_pts))

            if include_sample_times:
                xt_i_extra = data(xt)[groups_to_use[i], 0:data(ni)[groups_to_use[i]]]
                model_switch_i_extra = data(model_switch)[groups_to_use[i], 0:data(ni)[groups_to_use[i]]]
                # omitted
                tmp_include = ~(data(xt_i_extra) in data(xt_i))
                xt_i = np.concatenate(data(xt_i), data(xt_i_extra)[tmp_include])
                model_switch_i = np.concatenate(model_switch_i, model_switch_i_extra[tmp_include])
                tmp_order = np.argsort(xt_i)
                xt_i = xt_i[tmp_order]
                model_switch_i = model_switch_i[tmp_order]
        
        pred = None
        group_df = {"Time": xt_i, "PRED": pred, "Group": groups_to_use[i], "Model": model_switch_i}

        if predictions:
            bpop_val = parameters["bpop"][:, 1]
            b_ind = np.zeros([1, d_size])
            bocc_ind = np.zeros([docc_size, NumOcc])
            g0 = feval(model["fg_pointer"], x_i, a_i, bpop_val, b_ind, bocc_ind)

            pred_list = feval(model["ff_pointer"], model_switch_i, xt_i, g0, pypkpd_db)
            pred_list["pypkpd_db"] = None
            tmp = size(pred_list[0])
            filter(lambda a: a != 1, tmp)
            pred = np.array(pred_list[0]).reshape(tmp)
            group_df["PRED"] = pred
            
            if length(pred_list) > 1:
                extra_df = pred_list[-1]
                group_df["extra_df"] = extra_df

            # omitted
            # if PI:
                # sigma_full = get_dict_value(parameters, "sigma")
                # d_full = getfulld(get_dict_value(parameters, "docc")[:,1], get_dict_value(parameters, "covdocc"))

        if (include_a and a_i is not None):
            a_i = Matrix(a_i, axisnam=[None, None])
            group_df["a_i"] = a_i
        if (include_x and x_i is not None):
            x_i = Matrix(x_i, axisnam=[None, None])
            group_df["x_i"] = x_i
        
        # omitted IPRED
        df = group_df
    
    first_names = ["ID", "Time", "DV", "IPRED", "PRED"]
    tmp = []
    for i in first_names:
        if i in df.keys():
            tmp.append(i)
    first_names = tmp
    tmp = []
    for i in df.keys():
        if i not in first_names:
            tmp.append(i)
    other_names = tmp

    # omitted

    return df