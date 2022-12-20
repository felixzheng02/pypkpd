"""
## Create global variables in the poped database
## 
## Function takes design variables from input files
## and converts them to the global variables needed
## in poped.  Typically not used by the user.  Instead 
## use the function \code{\link{create.poped.database}}.
## 
## @param pypkpd_db A PopED database 
## @return A PopED database
## @family poped_input
## @export
## @keywords internal

## Function written to match MATLAB function convert_variables()


## Author: Caiya Zhang, Yuchen Zheng
"""


def convert_variables(pypkpd_db):
    design = pypkpd_db["design"]
    design_space = pypkpd_db["design_space"]

    #pypkpd_db$ga = zeros(pypkpd_db["design"]m,size(pypkpd_db["design"]a,2))
    #pypkpd_db$gx = zeros(pypkpd_db["design"]m,size(pypkpd_db["design"]x,2))
    #pypkpd_db$gxt = zeros(pypkpd_db["design"]m,max(pypkpd_db["design_space"]maxni))
    #pypkpd_db$gminxt = zeros(pypkpd_db["design"]m,max(pypkpd_db["design_space"]maxni))
    #pypkpd_db$gmaxxt = zeros(pypkpd_db["design"]m,max(pypkpd_db["design_space"]maxni))
    #pypkpd_db$gmaxa = zeros(pypkpd_db["design"]m,size(pypkpd_db["design"]a,2))
    #pypkpd_db$gmina = zeros(pypkpd_db["design"]m,size(pypkpd_db["design"]a,2))

    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, 1),nrow=1,byrow=T),design$ni,'ni')==1)){
    #         pypkpd_db$gni = design$ni
    #     }    
    
    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, max(pypkpd_db["design_space"]maxni)),nrow=1,byrow=T),design$model_switch,'model_switch')==1)){
    #         if(is.null(pypkpd_db$global_model_switch)){
    #             pypkpd_db$global_model_switch=design$model_switch
    #         } else {
    #             pypkpd_db$global_model_switch[1:pypkpd_db["design"]m,1:max(pypkpd_db["design_space"]maxni)]=design$model_switch
    #         }
    #     }

    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, max(pypkpd_db["design_space"]maxni)),nrow=1,byrow=T),design$xt,'xt')==1)){
    #         pypkpd_db$gxt[1:pypkpd_db["design"]m,1:max(pypkpd_db["design_space"]maxni)]=design$xt
    #     }
    
    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, max(pypkpd_db["design_space"]maxni)),nrow=1,byrow=T),design_space$maxxt,'maxxt')==1)){
    #         pypkpd_db$gmaxxt[1:pypkpd_db["design"]m,1:max(pypkpd_db["design_space"]maxni)] = design_space$maxxt
    #     }

    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, max(pypkpd_db["design_space"]maxni)),nrow=1,byrow=T),design_space$minxt,'minxt')==1)){
    #         pypkpd_db$gminxt[1:pypkpd_db["design"]m,1:max(pypkpd_db["design_space"]maxni)] = design_space$minxt
    #     }
    
    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, size(pypkpd_db["design"]a,2)),nrow=1,byrow=T),design$a,'a')==1)){
    #       pypkpd_db$ga[1:pypkpd_db["design"]m,1:size(pypkpd_db["design"]a,2)]=design$a
    #     }
       
    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, size(pypkpd_db["design"]a,2)),nrow=1,byrow=T),design_space$maxa,'maxa')==1)){
    #         pypkpd_db$gmaxa[1:pypkpd_db["design"]m,1:size(pypkpd_db["design"]a,2)]=design_space$maxa
    #     }
    
    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, size(pypkpd_db["design"]a,2)),nrow=1,byrow=T),design_space$mina,'mina')==1)){
    #         pypkpd_db$gmina[1:pypkpd_db["design"]m,1:size(pypkpd_db["design"]a,2)]=design_space$mina
    #     }
    
    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, size(pypkpd_db["design"]x,2)),nrow=1,byrow=T),design$x,'x')==1)){
    #         pypkpd_db$gx[1:pypkpd_db["design"]m,1:size(pypkpd_db["design"]x,2)]=design$x
    #     }
    
    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, 1),nrow=1,byrow=T),design$groupsize,'groupsize')==1)){
    #         pypkpd_db$groupsize=design$groupsize
    #     }
    
    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, max(pypkpd_db["design_space"]maxni)),nrow=1,byrow=T),
    #                       design_space$G_xt,'G')==1)){
    #         pypkpd_db$G=design_space$G_xt
    #     }
    
    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, size(pypkpd_db["design"]a,2)),nrow=1,byrow=T),design_space$G_a,'Ga')==1)){
    #         if(is.null(dim(design_space$G_a))){
    #             design_space$G_a=Matrix(data=design_space$G_a,nrow=1,byrow=True)
    #         }
    #         pypkpd_db$Ga=design_space$G_a
    #     }
    
    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, size(pypkpd_db["design"]x,2)),nrow=1,byrow=T),design_space$G_x,'Gx')==1)){
    #         pypkpd_db$Gx=design_space$G_x
    #     }
    
    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, 1),nrow=1,byrow=T),design_space$maxgroupsize,'maxgroupsize')==1)){
    #         pypkpd_db$maxgroupsize=design_space$maxgroupsize
    #     }

    #     if((test_mat_size(Matrix(c(pypkpd_db["design"]m, 1),nrow=1,byrow=T),design_space$mingroupsize,'mingroupsize')==1)){
    #         pypkpd_db$mingroupsize=design_space$mingroupsize
    #     }
    
    return pypkpd_db
