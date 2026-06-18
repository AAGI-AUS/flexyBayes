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
      [1] 54
      
      $.refusal_registry$keys
       [1] "approximate_route_not_yet_registered"         
       [2] "approximation_absent"                         
       [3] "approximation_no_smooth_path"                 
       [4] "approximation_scheme_unknown"                 
       [5] "approximation_spec_invalid"                   
       [6] "block_not_positive_definite"                  
       [7] "block_partition_incomplete"                   
       [8] "blocks_empty_list"                            
       [9] "blocks_not_a_list"                            
      [10] "blocks_not_in_known_matrices"                 
      [11] "chol_not_in_known_matrices"                   
      [12] "chol_not_square"                              
      [13] "chol_not_triangular"                          
      [14] "code_flags_mutually_exclusive"                
      [15] "cov_arg_not_fb_cov"                           
      [16] "design_memory_exceeds_ceiling"                
      [17] "engine_pin_backend_conflict"                  
      [18] "fa_rank_exceeds_dim"                          
      [19] "fa_rank_invalid"                              
      [20] "fb_cov_missing_matrix"                        
      [21] "fb_cov_type_unknown"                          
      [22] "formula_not_two_sided"                        
      [23] "grammar_brms_known_matrices_unsupported"      
      [24] "grammar_brms_with_asreml_terms"               
      [25] "gretaR_below_version_floor"                   
      [26] "gretaR_cannot_represent_structured_cov"       
      [27] "gretaR_family_unsupported"                    
      [28] "gretaR_not_installed"                         
      [29] "gretaR_random_group_not_in_data"              
      [30] "gretaR_random_term_type_unsupported"          
      [31] "heterogeneous_residual_factor_not_in_cell_key"
      [32] "known_matrices_data_name_collision"           
      [33] "known_matrix_dim_mismatch"                    
      [34] "known_matrix_dimnames_mismatch"               
      [35] "known_matrix_level_mismatch"                  
      [36] "low_rank_rank_exceeds_basis"                  
      [37] "low_rank_rank_invalid"                        
      [38] "low_rank_requires_greta"                      
      [39] "low_rank_scheme_required"                     
      [40] "native_greta_requires_greta_backend"          
      [41] "precision_not_in_known_matrices"              
      [42] "precision_not_positive_definite"              
      [43] "precision_not_square"                         
      [44] "precision_not_symmetric"                      
      [45] "predict_kernel_invalid_include"               
      [46] "rcov_type_unsupported_for_aggregation"        
      [47] "representation_unknown_for_preflight"         
      [48] "response_not_in_data"                         
      [49] "review_code_backend_unsupported"              
      [50] "smooth_variable_not_in_data"                  
      [51] "stan_cannot_represent_structured_cov"         
      [52] "tensor_smooth_unsupported"                    
      [53] "unsupported_family"                           
      [54] "vm_redundant_specification"                   
      
      
      $.backend_registry
      $.backend_registry$registry
      [1] ".backend_registry"
      
      $.backend_registry$locked
      [1] TRUE
      
      $.backend_registry$n_keys
      [1] 4
      
      $.backend_registry$keys
      [1] "brms"   "greta"  "gretaR" "inla"  
      
      

