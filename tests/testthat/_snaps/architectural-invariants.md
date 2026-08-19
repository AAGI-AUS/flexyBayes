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
      [1] 116
      
      $.refusal_registry$keys
        [1] "aggregation_response_incomplete"              
        [2] "aggregation_route_unavailable"                
        [3] "approximate_route_not_yet_registered"         
        [4] "approximation_absent"                         
        [5] "approximation_no_smooth_path"                 
        [6] "approximation_scheme_unknown"                 
        [7] "approximation_spec_invalid"                   
        [8] "ar1_residual_not_representable"               
        [9] "ar1_spatial_requires_complete_grid"           
       [10] "ar2_not_representable"                        
       [11] "asreml_function_not_recognised"               
       [12] "at_level_conditioning_unsupported"            
       [13] "augment_cell_not_determinable"                
       [14] "auto_no_active_route"                         
       [15] "backend_quarantined"                          
       [16] "block_not_positive_definite"                  
       [17] "block_partition_incomplete"                   
       [18] "blocks_empty_list"                            
       [19] "blocks_not_a_list"                            
       [20] "blocks_not_in_known_matrices"                 
       [21] "brms_cannot_augment_nongaussian"              
       [22] "brms_cannot_represent_term"                   
       [23] "brms_factor_random_slope_unsupported"         
       [24] "brms_ingest_feature_unsupported"              
       [25] "brms_random_effect_form_unsupported"          
       [26] "chol_not_in_known_matrices"                   
       [27] "chol_not_square"                              
       [28] "chol_not_triangular"                          
       [29] "classify_random_factor_not_supported"         
       [30] "classify_requires_emmeans"                    
       [31] "code_flags_mutually_exclusive"                
       [32] "conditional_loglik_not_available"             
       [33] "corh_no_equicorrelation_representation"       
       [34] "cov_arg_not_fb_cov"                           
       [35] "covariate_zero_fill_not_supported"            
       [36] "design_memory_exceeds_ceiling"                
       [37] "dsum_structured_inner_unsupported"            
       [38] "engine_pin_backend_conflict"                  
       [39] "fa_not_representable"                         
       [40] "fa_rank_exceeds_dim"                          
       [41] "fa_rank_invalid"                              
       [42] "family_argument_not_recognised"               
       [43] "fb_cov_missing_matrix"                        
       [44] "fb_cov_type_unknown"                          
       [45] "fit_lacks_posterior_draws"                    
       [46] "fixed_smoother_not_supported"                 
       [47] "formula_not_two_sided"                        
       [48] "grammar_brms_known_matrices_unsupported"      
       [49] "grammar_brms_with_asreml_terms"               
       [50] "gretaR_below_version_floor"                   
       [51] "gretaR_cannot_represent_structured_cov"       
       [52] "gretaR_family_unsupported"                    
       [53] "gretaR_not_installed"                         
       [54] "gretaR_random_group_not_in_data"              
       [55] "gretaR_random_term_type_unsupported"          
       [56] "heterogeneous_residual_factor_not_in_cell_key"
       [57] "heterogeneous_residual_family_has_no_sigma"   
       [58] "heterogeneous_residual_multiple_factors"      
       [59] "inla_gate_refused"                            
       [60] "inla_variable_used_twice"                     
       [61] "interaction_not_representable"                
       [62] "known_matrices_data_name_collision"           
       [63] "known_matrix_dim_mismatch"                    
       [64] "known_matrix_dimnames_mismatch"               
       [65] "known_matrix_level_mismatch"                  
       [66] "loo_requires_sampler_draws"                   
       [67] "low_rank_rank_exceeds_basis"                  
       [68] "low_rank_rank_invalid"                        
       [69] "low_rank_requires_greta"                      
       [70] "low_rank_scheme_required"                     
       [71] "met_summary_not_available"                    
       [72] "missing_covariate_not_supported"              
       [73] "missing_response_refused"                     
       [74] "model_matrix_not_recoverable"                 
       [75] "native_greta_fit_quarantined"                 
       [76] "native_greta_requires_greta_backend"          
       [77] "numeric_variable_in_random_interaction"       
       [78] "pp_check_requires_predictive_draws"           
       [79] "precision_not_in_known_matrices"              
       [80] "precision_not_positive_definite"              
       [81] "precision_not_square"                         
       [82] "precision_not_symmetric"                      
       [83] "predict_kernel_invalid_include"               
       [84] "prior_argument_duplicated"                    
       [85] "prior_argument_missing"                       
       [86] "prior_argument_unknown"                       
       [87] "prior_distribution_not_a_call"                
       [88] "prior_distribution_unknown"                   
       [89] "prior_hyperparameter_not_scalar"              
       [90] "prior_hyperparameter_out_of_domain"           
       [91] "prior_not_translatable_for_backend"           
       [92] "prior_spec_empty"                             
       [93] "prior_spec_not_formula"                       
       [94] "prior_spec_not_two_sided"                     
       [95] "prior_target_argument_missing"                
       [96] "prior_target_not_in_model"                    
       [97] "prior_target_unsupported"                     
       [98] "representation_unknown_for_preflight"         
       [99] "residual_type_unsupported_for_aggregation"    
      [100] "response_not_in_data"                         
      [101] "review_code_backend_unsupported"              
      [102] "row_count_exceeds_integer"                    
      [103] "smooth_variable_not_in_data"                  
      [104] "stan_cannot_represent_ar1_field"              
      [105] "stan_cannot_represent_ar1_residual"           
      [106] "stan_cannot_represent_structured_cov"         
      [107] "stan_cannot_represent_structured_residual"    
      [108] "str_not_representable"                        
      [109] "tensor_smooth_unsupported"                    
      [110] "term_in_fixed_and_random"                     
      [111] "triangulate_incomparable_fits"                
      [112] "unsupported_family"                           
      [113] "update_call_not_reconstructable"              
      [114] "variogram_requires_design_index"              
      [115] "vm_redundant_specification"                   
      [116] "weights_not_supported"                        
      
      
      $.backend_registry
      $.backend_registry$registry
      [1] ".backend_registry"
      
      $.backend_registry$locked
      [1] TRUE
      
      $.backend_registry$n_keys
      [1] 4
      
      $.backend_registry$keys
      [1] "brms"   "greta"  "gretaR" "inla"  
      
      

