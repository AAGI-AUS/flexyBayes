# the verb list is stable (snapshot)

    Code
      sort(.FB_VERBS)
    Output
      [1] "fb_brms"                "fb_inla"                "fb_plan"               
      [4] "flexybayes"             "validate_approximation"

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
      [1] 123
      
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
       [16] "binomial_response_not_binary"                 
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
       [27] "brms_stan_toolchain_compile_failure"          
       [28] "cell_count_exceeds_integer"                   
       [29] "chol_not_in_known_matrices"                   
       [30] "chol_not_square"                              
       [31] "chol_not_triangular"                          
       [32] "classify_random_factor_not_supported"         
       [33] "classify_requires_emmeans"                    
       [34] "code_flags_mutually_exclusive"                
       [35] "conditional_loglik_not_available"             
       [36] "cor_slope_prior_not_translatable_inla"        
       [37] "corh_no_equicorrelation_representation"       
       [38] "cov_arg_not_fb_cov"                           
       [39] "covariate_zero_fill_not_supported"            
       [40] "design_memory_exceeds_ceiling"                
       [41] "dsum_structured_inner_unsupported"            
       [42] "engine_pin_backend_conflict"                  
       [43] "fa_not_representable"                         
       [44] "fa_rank_exceeds_dim"                          
       [45] "fa_rank_invalid"                              
       [46] "family_argument_not_recognised"               
       [47] "fb_cov_missing_matrix"                        
       [48] "fb_cov_type_unknown"                          
       [49] "fit_lacks_posterior_draws"                    
       [50] "fixed_smoother_not_supported"                 
       [51] "formula_not_two_sided"                        
       [52] "grammar_brms_known_matrices_unsupported"      
       [53] "grammar_brms_with_asreml_terms"               
       [54] "heterogeneous_residual_factor_not_in_cell_key"
       [55] "heterogeneous_residual_family_has_no_sigma"   
       [56] "heterogeneous_residual_multiple_factors"      
       [57] "inla_gate_refused"                            
       [58] "inla_program_failed"                          
       [59] "inla_variable_used_twice"                     
       [60] "interaction_not_representable"                
       [61] "known_matrices_data_name_collision"           
       [62] "known_matrix_dim_mismatch"                    
       [63] "known_matrix_dimnames_mismatch"               
       [64] "known_matrix_level_mismatch"                  
       [65] "loo_requires_sampler_draws"                   
       [66] "low_rank_rank_exceeds_basis"                  
       [67] "low_rank_rank_invalid"                        
       [68] "low_rank_scheme_required"                     
       [69] "low_rank_smooth_unsupported"                  
       [70] "met_summary_not_available"                    
       [71] "missing_covariate_not_supported"              
       [72] "missing_response_refused"                     
       [73] "model_matrix_not_recoverable"                 
       [74] "nested_index_name_collision"                  
       [75] "nested_variable_not_in_data"                  
       [76] "numeric_variable_in_random_interaction"       
       [77] "pp_check_requires_predictive_draws"           
       [78] "precision_not_in_known_matrices"              
       [79] "precision_not_positive_definite"              
       [80] "precision_not_square"                         
       [81] "precision_not_symmetric"                      
       [82] "prior_argument_duplicated"                    
       [83] "prior_argument_missing"                       
       [84] "prior_argument_unknown"                       
       [85] "prior_distribution_not_a_call"                
       [86] "prior_distribution_unknown"                   
       [87] "prior_hyperparameter_not_scalar"              
       [88] "prior_hyperparameter_out_of_domain"           
       [89] "prior_not_translatable_for_backend"           
       [90] "prior_spec_empty"                             
       [91] "prior_spec_not_formula"                       
       [92] "prior_spec_not_two_sided"                     
       [93] "prior_target_argument_missing"                
       [94] "prior_target_not_in_model"                    
       [95] "prior_target_unsupported"                     
       [96] "representation_unknown_for_preflight"         
       [97] "residual_type_unsupported_for_aggregation"    
       [98] "response_not_in_data"                         
       [99] "review_code_backend_unsupported"              
      [100] "row_count_exceeds_integer"                    
      [101] "smooth_variable_not_in_data"                  
      [102] "stan_cannot_represent_ar1_field"              
      [103] "stan_cannot_represent_ar1_residual"           
      [104] "stan_cannot_represent_structured_cov"         
      [105] "stan_cannot_represent_structured_residual"    
      [106] "str_not_representable"                        
      [107] "stream_response_must_be_a_column"             
      [108] "tensor_smooth_unsupported"                    
      [109] "term_in_fixed_and_random"                     
      [110] "triangulate_incomparable_fits"                
      [111] "unknown_backend"                              
      [112] "unsupported_family"                           
      [113] "update_call_not_reconstructable"              
      [114] "update_unnamed_argument_not_supported"        
      [115] "us_gxe_prior_not_translatable_inla"           
      [116] "variogram_requires_design_index"              
      [117] "vm_redundant_specification"                   
      [118] "vpredict_component_not_in_posterior"          
      [119] "vpredict_functional_not_one_value_per_draw"   
      [120] "vpredict_no_functional"                       
      [121] "weights_not_aggregatable"                     
      [122] "weights_not_supported"                        
      [123] "weights_requires_gaussian"                    
      
      
      $.backend_registry
      $.backend_registry$registry
      [1] ".backend_registry"
      
      $.backend_registry$locked
      [1] TRUE
      
      $.backend_registry$n_keys
      [1] 2
      
      $.backend_registry$keys
      [1] "brms" "inla"
      
      

