"""
Step-by-step PyDESeq2 workflow
===============================

This notebook details all the steps of the PyDESeq2 pipeline.

It allows you to run the PyDESeq2 pipeline on the synthetic data provided
in this repository.

If this is your first contact with PyDESeq2, we recommend you first have a look at the
:doc:`standard workflow example <plot_minimal_pydeseq2_pipeline>`.

.. contents:: Contents
    :local:
    :depth: 3

We start by importing required packages and setting up an optional path to save
results.
"""

import os
import pickle as pkl

import numpy as np

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data

save = False                                                

if save:
                                                                                
           
    output_path = "../output_files/synthetic_example"
    os.makedirs(output_path, exist_ok=True)                                   

    
              
              
 
                                                                              
                                                                                   
         
 
                                                                                 
                                          
                                                                
                                                                            
                                   
 
                                               
                                                                         
 
                             
 
                                                                             
             
                                                                        
                                           

counts_df = load_example_data(
    modality="raw_counts",
    dataset="synthetic",
    debug=False,
)

metadata = load_example_data(
    modality="metadata",
    dataset="synthetic",
    debug=False,
)

    
                         
                         
                                                    
                           
 
                                                                      
                           
                                                                            
 
                                       
                                                                                   
                                                          
                                                                                   
 
           
                                                                        
             

inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design="~condition",                                            
                         
    refit_cooks=True,
    inference=inference,
)

    
                               
                               

dds.fit_size_factors()

dds.obs["size_factors"]

    
                          
                          

dds.fit_genewise_dispersions()

dds.var["genewise_dispersions"]

    
                                   
                                   

dds.fit_dispersion_trend()
dds.uns["trend_coeffs"]
dds.var["fitted_dispersions"]

    
                   
                   

dds.fit_dispersion_prior()
print(
    f"logres_prior={dds.uns['_squared_logres']}, sigma_prior={dds.uns['prior_disp_var']}"
)

    
                 
                 
                                                                         
                       
                                                                             
                                 
                                                                              
                              

dds.fit_MAP_dispersions()
dds.var["MAP_dispersions"]
dds.var["dispersions"]

    
                      
                      
                                                                            
         
                                                                             
                                                          

dds.fit_LFC()
dds.varm["LFC"]

    
                                     
                                     
                                  

dds.calculate_cooks()
if dds.refit_cooks:
                            
    dds.refit()

    
                 

if save:
    with open(os.path.join(output_path, "dds_detailed_pipe.pkl"), "wb") as f:
        pkl.dump(dds, f)

    
                         
                         
                                                                          
                                                                                  
                                                                               
                                                                                 
                                                                 
                         
 
                                                                    
                                                                      
                                                                             
                   

ds = DeseqStats(
    dds,
    contrast=np.array([0, 1]),
    alpha=0.05,
    cooks_filter=True,
    independent_filter=True,
)

    
            
            

ds.run_wald_test()
ds.p_values

    
                 
                 
                  
 
           
                                                                        
             

if ds.cooks_filter:
    ds._cooks_filtering()
ds.p_values

    
                    
                    

if ds.independent_filter:
    ds._independent_filtering()
else:
    ds._p_value_adjustment()

ds.padj

    
                              
                              
                                                                            
        

ds.summary()

    
                                        

if save:
    with open(os.path.join(output_path, "stat_results_detailed_pipe.pkl"), "wb") as f:
        pkl.dump(ds, f)

    
               
               
                                                                                
                                                                

ds.lfc_shrink(coeff="condition[T.B]")


    
                 

if save:
    with open(os.path.join(output_path, "shrunk_results_detailed_pipe.pkl"), "wb") as f:
        pkl.dump(ds, f)
