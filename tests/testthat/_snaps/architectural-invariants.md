# the verb list is stable (snapshot)

    Code
      sort(.FB_VERBS)
    Output
      [1] "fb_brms"                "fb_greta"               "fb_inla"               
      [4] "fb_plan"                "flexybayes"             "triangulate"           
      [7] "validate_approximation"

# the constructor schemas are stable (snapshot)

    Code
      lapply(inst, .fb_schema_dump)
    Output
      $fb_prior
      $fb_prior$class
      [1] "fb_prior" "list"    
      
      $fb_prior$elements
      [1] "specs"
      
      $fb_prior$attributes
      character(0)
      
      
      $fb_cov
      $fb_cov$class
      [1] "fb_cov" "list"  
      
      $fb_cov$elements
      [1] "M"      "levels" "scheme" "type"  
      
      $fb_cov$attributes
      [1] "representation_class" "type"                 "validation_summary"  
      
      
      $fb_approx
      $fb_approx$class
      [1] "fb_approx" "list"     
      
      $fb_approx$elements
      [1] "rank"   "scheme"
      
      $fb_approx$attributes
      [1] "bias_bound_promise"
      
      
      $fb_engine
      $fb_engine$class
      [1] "fb_engine" "list"     
      
      $fb_engine$elements
      [1] "name"             "opts"             "paradigm"         "toolchain_status"
      
      $fb_engine$attributes
      character(0)
      
      

# the registry key inventory is stable (snapshot)

    Code
      dump
    Output
      $.representation_registry
      $.representation_registry$registry
      [1] ".representation_registry"
      
      $.representation_registry$locked
      [1] TRUE
      
      $.representation_registry$n_keys
      [1] 16
      
      $.representation_registry$keys
       [1] "banded_smooth"                "block_diagonal"              
       [3] "chol_cov"                     "dense_baseline"              
       [5] "dense_cov"                    "dense_smooth"                
       [7] "indexed_fixed_factor"         "indexed_fixed_factor_numeric"
       [9] "indexed_fixed_numeric"        "indexed_random_intercept"    
      [11] "indexed_structured_estimate"  "indexed_structured_known"    
      [13] "low_rank"                     "pedigree_sparse_precision"   
      [15] "sparse_precision"             "sparse_smooth"               
      
      
      $.approximation_registry
      $.approximation_registry$registry
      [1] ".approximation_registry"
      
      $.approximation_registry$locked
      [1] TRUE
      
      $.approximation_registry$n_keys
      [1] 1
      
      $.approximation_registry$keys
      [1] "low_rank_smooth"
      
      
      $.backend_independence_registry
      $.backend_independence_registry$registry
      [1] ".backend_independence_registry"
      
      $.backend_independence_registry$locked
      [1] TRUE
      
      $.backend_independence_registry$n_keys
      [1] 3
      
      $.backend_independence_registry$keys
      [1] "brms||greta" "brms||inla"  "greta||inla"
      
      
      $.refusal_registry
      $.refusal_registry$registry
      [1] ".refusal_registry"
      
      $.refusal_registry$locked
      [1] TRUE
      
      $.refusal_registry$n_keys
      [1] 88
      
      $.refusal_registry$keys
       [1] "approximate_route_not_yet_registered"         
       [2] "approximation_absent"                         
       [3] "approximation_no_smooth_path"                 
       [4] "approximation_scheme_unknown"                 
       [5] "approximation_spec_invalid"                   
       [6] "ar1_residual_not_representable"               
       [7] "ar1_spatial_requires_complete_grid"           
       [8] "ar2_not_representable"                        
       [9] "asreml_function_not_recognised"               
      [10] "at_level_conditioning_unsupported"            
      [11] "augment_cell_not_determinable"                
      [12] "auto_no_active_route"                         
      [13] "backend_quarantined"                          
      [14] "block_not_positive_definite"                  
      [15] "block_partition_incomplete"                   
      [16] "blocks_empty_list"                            
      [17] "blocks_not_a_list"                            
      [18] "blocks_not_in_known_matrices"                 
      [19] "brms_cannot_augment_nongaussian"              
      [20] "brms_cannot_represent_term"                   
      [21] "chol_not_in_known_matrices"                   
      [22] "chol_not_square"                              
      [23] "chol_not_triangular"                          
      [24] "code_flags_mutually_exclusive"                
      [25] "conditional_loglik_not_available"             
      [26] "corh_no_equicorrelation_representation"       
      [27] "cov_arg_not_fb_cov"                           
      [28] "design_memory_exceeds_ceiling"                
      [29] "dsum_structured_inner_unsupported"            
      [30] "engine_pin_backend_conflict"                  
      [31] "fa_not_representable"                         
      [32] "fa_rank_exceeds_dim"                          
      [33] "fa_rank_invalid"                              
      [34] "family_argument_not_recognised"               
      [35] "fb_cov_missing_matrix"                        
      [36] "fb_cov_type_unknown"                          
      [37] "fit_lacks_posterior_draws"                    
      [38] "formula_not_two_sided"                        
      [39] "grammar_brms_known_matrices_unsupported"      
      [40] "grammar_brms_with_asreml_terms"               
      [41] "gretaR_below_version_floor"                   
      [42] "gretaR_cannot_represent_structured_cov"       
      [43] "gretaR_family_unsupported"                    
      [44] "gretaR_not_installed"                         
      [45] "gretaR_random_group_not_in_data"              
      [46] "gretaR_random_term_type_unsupported"          
      [47] "heterogeneous_residual_factor_not_in_cell_key"
      [48] "heterogeneous_residual_family_has_no_sigma"   
      [49] "heterogeneous_residual_multiple_factors"      
      [50] "inla_gate_refused"                            
      [51] "interaction_not_representable"                
      [52] "known_matrices_data_name_collision"           
      [53] "known_matrix_dim_mismatch"                    
      [54] "known_matrix_dimnames_mismatch"               
      [55] "known_matrix_level_mismatch"                  
      [56] "low_rank_rank_exceeds_basis"                  
      [57] "low_rank_rank_invalid"                        
      [58] "low_rank_requires_greta"                      
      [59] "low_rank_scheme_required"                     
      [60] "met_summary_not_available"                    
      [61] "missing_covariate_not_supported"              
      [62] "missing_response_refused"                     
      [63] "model_matrix_not_recoverable"                 
      [64] "native_greta_fit_quarantined"                 
      [65] "native_greta_requires_greta_backend"          
      [66] "numeric_variable_in_random_interaction"       
      [67] "precision_not_in_known_matrices"              
      [68] "precision_not_positive_definite"              
      [69] "precision_not_square"                         
      [70] "precision_not_symmetric"                      
      [71] "predict_kernel_invalid_include"               
      [72] "representation_unknown_for_preflight"         
      [73] "residual_type_unsupported_for_aggregation"    
      [74] "response_not_in_data"                         
      [75] "review_code_backend_unsupported"              
      [76] "row_count_exceeds_integer"                    
      [77] "smooth_variable_not_in_data"                  
      [78] "stan_cannot_represent_ar1_field"              
      [79] "stan_cannot_represent_ar1_residual"           
      [80] "stan_cannot_represent_structured_cov"         
      [81] "stan_cannot_represent_structured_residual"    
      [82] "str_not_representable"                        
      [83] "tensor_smooth_unsupported"                    
      [84] "triangulate_incomparable_fits"                
      [85] "unsupported_family"                           
      [86] "update_call_not_reconstructable"              
      [87] "vm_redundant_specification"                   
      [88] "weights_not_supported"                        
      
      
      $.backend_registry
      $.backend_registry$registry
      [1] ".backend_registry"
      
      $.backend_registry$locked
      [1] TRUE
      
      $.backend_registry$n_keys
      [1] 4
      
      $.backend_registry$keys
      [1] "brms"   "greta"  "gretaR" "inla"  
      
      

