# the verb list is stable (snapshot)

    Code
      sort(.FB_VERBS)
    Output
      [1] "fb_brms"                "fb_inla"                "fb_plan"               
      [4] "flexybayes"             "triangulate"            "validate_approximation"

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
      [1] 1
      
      $.backend_independence_registry$keys
      [1] "brms||inla"
      
      
      $.refusal_registry
      $.refusal_registry$registry
      [1] ".refusal_registry"
      
      $.refusal_registry$locked
      [1] TRUE
      
      $.refusal_registry$n_keys
      [1] 112
      
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
       [26] "cell_count_exceeds_integer"                   
       [27] "chol_not_in_known_matrices"                   
       [28] "chol_not_square"                              
       [29] "chol_not_triangular"                          
       [30] "classify_random_factor_not_supported"         
       [31] "classify_requires_emmeans"                    
       [32] "code_flags_mutually_exclusive"                
       [33] "conditional_loglik_not_available"             
       [34] "corh_no_equicorrelation_representation"       
       [35] "cov_arg_not_fb_cov"                           
       [36] "covariate_zero_fill_not_supported"            
       [37] "design_memory_exceeds_ceiling"                
       [38] "dsum_structured_inner_unsupported"            
       [39] "engine_pin_backend_conflict"                  
       [40] "fa_not_representable"                         
       [41] "fa_rank_exceeds_dim"                          
       [42] "fa_rank_invalid"                              
       [43] "family_argument_not_recognised"               
       [44] "fb_cov_missing_matrix"                        
       [45] "fb_cov_type_unknown"                          
       [46] "fit_lacks_posterior_draws"                    
       [47] "fixed_smoother_not_supported"                 
       [48] "formula_not_two_sided"                        
       [49] "grammar_brms_known_matrices_unsupported"      
       [50] "grammar_brms_with_asreml_terms"               
       [51] "heterogeneous_residual_factor_not_in_cell_key"
       [52] "heterogeneous_residual_family_has_no_sigma"   
       [53] "heterogeneous_residual_multiple_factors"      
       [54] "inla_gate_refused"                            
       [55] "inla_program_failed"                          
       [56] "inla_variable_used_twice"                     
       [57] "interaction_not_representable"                
       [58] "known_matrices_data_name_collision"           
       [59] "known_matrix_dim_mismatch"                    
       [60] "known_matrix_dimnames_mismatch"               
       [61] "known_matrix_level_mismatch"                  
       [62] "loo_requires_sampler_draws"                   
       [63] "low_rank_rank_exceeds_basis"                  
       [64] "low_rank_rank_invalid"                        
       [65] "low_rank_scheme_required"                     
       [66] "low_rank_smooth_unsupported"                  
       [67] "met_summary_not_available"                    
       [68] "missing_covariate_not_supported"              
       [69] "missing_response_refused"                     
       [70] "model_matrix_not_recoverable"                 
       [71] "numeric_variable_in_random_interaction"       
       [72] "pp_check_requires_predictive_draws"           
       [73] "precision_not_in_known_matrices"              
       [74] "precision_not_positive_definite"              
       [75] "precision_not_square"                         
       [76] "precision_not_symmetric"                      
       [77] "prior_argument_duplicated"                    
       [78] "prior_argument_missing"                       
       [79] "prior_argument_unknown"                       
       [80] "prior_distribution_not_a_call"                
       [81] "prior_distribution_unknown"                   
       [82] "prior_hyperparameter_not_scalar"              
       [83] "prior_hyperparameter_out_of_domain"           
       [84] "prior_not_translatable_for_backend"           
       [85] "prior_spec_empty"                             
       [86] "prior_spec_not_formula"                       
       [87] "prior_spec_not_two_sided"                     
       [88] "prior_target_argument_missing"                
       [89] "prior_target_not_in_model"                    
       [90] "prior_target_unsupported"                     
       [91] "representation_unknown_for_preflight"         
       [92] "residual_type_unsupported_for_aggregation"    
       [93] "response_not_in_data"                         
       [94] "review_code_backend_unsupported"              
       [95] "row_count_exceeds_integer"                    
       [96] "smooth_variable_not_in_data"                  
       [97] "stan_cannot_represent_ar1_field"              
       [98] "stan_cannot_represent_ar1_residual"           
       [99] "stan_cannot_represent_structured_cov"         
      [100] "stan_cannot_represent_structured_residual"    
      [101] "str_not_representable"                        
      [102] "tensor_smooth_unsupported"                    
      [103] "term_in_fixed_and_random"                     
      [104] "triangulate_incomparable_fits"                
      [105] "unknown_backend"                              
      [106] "unsupported_family"                           
      [107] "update_call_not_reconstructable"              
      [108] "variogram_requires_design_index"              
      [109] "vm_redundant_specification"                   
      [110] "weights_not_aggregatable"                     
      [111] "weights_not_supported"                        
      [112] "weights_requires_gaussian"                    
      
      
      $.backend_registry
      $.backend_registry$registry
      [1] ".backend_registry"
      
      $.backend_registry$locked
      [1] TRUE
      
      $.backend_registry$n_keys
      [1] 2
      
      $.backend_registry$keys
      [1] "brms" "inla"
      
      

