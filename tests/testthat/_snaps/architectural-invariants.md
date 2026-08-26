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
      [1] 121
      
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
       [12] "at_field_per_level_hyper_not_representable"   
       [13] "at_level_conditioning_unsupported"            
       [14] "augment_cell_not_determinable"                
       [15] "auto_no_active_route"                         
       [16] "backend_quarantined"                          
       [17] "block_not_positive_definite"                  
       [18] "block_partition_incomplete"                   
       [19] "blocks_empty_list"                            
       [20] "blocks_not_a_list"                            
       [21] "blocks_not_in_known_matrices"                 
       [22] "brms_cannot_augment_nongaussian"              
       [23] "brms_cannot_represent_term"                   
       [24] "brms_factor_random_slope_unsupported"         
       [25] "brms_ingest_feature_unsupported"              
       [26] "brms_random_effect_form_unsupported"          
       [27] "cell_count_exceeds_integer"                   
       [28] "chol_not_in_known_matrices"                   
       [29] "chol_not_square"                              
       [30] "chol_not_triangular"                          
       [31] "classify_random_factor_not_supported"         
       [32] "classify_requires_emmeans"                    
       [33] "code_flags_mutually_exclusive"                
       [34] "conditional_loglik_not_available"             
       [35] "corh_no_equicorrelation_representation"       
       [36] "cov_arg_not_fb_cov"                           
       [37] "covariate_zero_fill_not_supported"            
       [38] "design_memory_exceeds_ceiling"                
       [39] "dsum_structured_inner_unsupported"            
       [40] "engine_pin_backend_conflict"                  
       [41] "fa_not_representable"                         
       [42] "fa_rank_exceeds_dim"                          
       [43] "fa_rank_invalid"                              
       [44] "family_argument_not_recognised"               
       [45] "fb_cov_missing_matrix"                        
       [46] "fb_cov_type_unknown"                          
       [47] "fit_lacks_posterior_draws"                    
       [48] "fixed_smoother_not_supported"                 
       [49] "formula_not_two_sided"                        
       [50] "grammar_brms_known_matrices_unsupported"      
       [51] "grammar_brms_with_asreml_terms"               
       [52] "gretaR_below_version_floor"                   
       [53] "gretaR_cannot_represent_structured_cov"       
       [54] "gretaR_family_unsupported"                    
       [55] "gretaR_not_installed"                         
       [56] "gretaR_random_group_not_in_data"              
       [57] "gretaR_random_term_type_unsupported"          
       [58] "heterogeneous_residual_factor_not_in_cell_key"
       [59] "heterogeneous_residual_family_has_no_sigma"   
       [60] "heterogeneous_residual_multiple_factors"      
       [61] "inla_gate_refused"                            
       [62] "inla_program_failed"                          
       [63] "inla_variable_used_twice"                     
       [64] "interaction_not_representable"                
       [65] "known_matrices_data_name_collision"           
       [66] "known_matrix_dim_mismatch"                    
       [67] "known_matrix_dimnames_mismatch"               
       [68] "known_matrix_level_mismatch"                  
       [69] "loo_requires_sampler_draws"                   
       [70] "low_rank_rank_exceeds_basis"                  
       [71] "low_rank_rank_invalid"                        
       [72] "low_rank_requires_greta"                      
       [73] "low_rank_scheme_required"                     
       [74] "met_summary_not_available"                    
       [75] "missing_covariate_not_supported"              
       [76] "missing_response_refused"                     
       [77] "model_matrix_not_recoverable"                 
       [78] "native_greta_fit_quarantined"                 
       [79] "native_greta_requires_greta_backend"          
       [80] "numeric_variable_in_random_interaction"       
       [81] "pp_check_requires_predictive_draws"           
       [82] "precision_not_in_known_matrices"              
       [83] "precision_not_positive_definite"              
       [84] "precision_not_square"                         
       [85] "precision_not_symmetric"                      
       [86] "predict_kernel_invalid_include"               
       [87] "prior_argument_duplicated"                    
       [88] "prior_argument_missing"                       
       [89] "prior_argument_unknown"                       
       [90] "prior_distribution_not_a_call"                
       [91] "prior_distribution_unknown"                   
       [92] "prior_hyperparameter_not_scalar"              
       [93] "prior_hyperparameter_out_of_domain"           
       [94] "prior_not_translatable_for_backend"           
       [95] "prior_spec_empty"                             
       [96] "prior_spec_not_formula"                       
       [97] "prior_spec_not_two_sided"                     
       [98] "prior_target_argument_missing"                
       [99] "prior_target_not_in_model"                    
      [100] "prior_target_unsupported"                     
      [101] "representation_unknown_for_preflight"         
      [102] "residual_type_unsupported_for_aggregation"    
      [103] "response_not_in_data"                         
      [104] "review_code_backend_unsupported"              
      [105] "row_count_exceeds_integer"                    
      [106] "smooth_variable_not_in_data"                  
      [107] "stan_cannot_represent_ar1_field"              
      [108] "stan_cannot_represent_ar1_residual"           
      [109] "stan_cannot_represent_structured_cov"         
      [110] "stan_cannot_represent_structured_residual"    
      [111] "str_not_representable"                        
      [112] "tensor_smooth_unsupported"                    
      [113] "term_in_fixed_and_random"                     
      [114] "triangulate_incomparable_fits"                
      [115] "unsupported_family"                           
      [116] "update_call_not_reconstructable"              
      [117] "variogram_requires_design_index"              
      [118] "vm_redundant_specification"                   
      [119] "weights_not_aggregatable"                     
      [120] "weights_not_supported"                        
      [121] "weights_requires_gaussian"                    
      
      
      $.backend_registry
      $.backend_registry$registry
      [1] ".backend_registry"
      
      $.backend_registry$locked
      [1] TRUE
      
      $.backend_registry$n_keys
      [1] 4
      
      $.backend_registry$keys
      [1] "brms"   "greta"  "gretaR" "inla"  
      
      

