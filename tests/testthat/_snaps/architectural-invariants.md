# the verb list is stable (snapshot)

    Code
      cat(sort(.FB_VERBS), sep = "\n")
    Output
      fb_brms
      fb_inla
      fb_plan
      flexybayes
      validate_approximation

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
      for (d in dump) {
        cat(d$registry, ": locked=", d$locked, " n_keys=", d$n_keys, "\n", sep = "")
        cat(paste0("  ", d$keys), sep = "\n")
        cat("\n")
      }
    Output
      .representation_registry: locked=TRUE n_keys=16
        banded_smooth
        block_diagonal
        chol_cov
        dense_baseline
        dense_cov
        dense_smooth
        indexed_fixed_factor
        indexed_fixed_factor_numeric
        indexed_fixed_numeric
        indexed_random_intercept
        indexed_structured_estimate
        indexed_structured_known
        low_rank
        pedigree_sparse_precision
        sparse_precision
        sparse_smooth
      
      .approximation_registry: locked=TRUE n_keys=1
        low_rank_smooth
      
      .backend_independence_registry: locked=TRUE n_keys=1
        brms||inla
      
      .refusal_registry: locked=TRUE n_keys=190
        aggregation_response_incomplete
        aggregation_route_unavailable
        approximate_route_not_yet_registered
        approximation_absent
        approximation_no_smooth_path
        approximation_scheme_unknown
        approximation_spec_invalid
        ar1_residual_not_representable
        ar1_spatial_requires_complete_grid
        ar2_not_representable
        asreml_function_not_recognised
        at_field_per_level_hyper_not_representable
        at_level_conditioning_unsupported
        augment_cell_not_determinable
        auto_no_active_route
        binomial_response_not_binary
        block_not_positive_definite
        block_partition_incomplete
        blocks_empty_list
        blocks_not_a_list
        blocks_not_in_known_matrices
        brms_cannot_augment_nongaussian
        brms_cannot_represent_term
        brms_factor_random_slope_unsupported
        brms_ingest_feature_unsupported
        brms_random_effect_form_unsupported
        brms_stan_toolchain_compile_failure
        cell_count_exceeds_integer
        chol_not_in_known_matrices
        chol_not_square
        chol_not_triangular
        classify_random_factor_not_supported
        classify_requires_emmeans
        code_flags_mutually_exclusive
        conditional_loglik_not_available
        cor_slope_prior_not_translatable_inla
        corh_no_equicorrelation_representation
        correlated_covariance_requires_brms
        cov_arg_not_fb_cov
        covariate_zero_fill_not_supported
        design_memory_exceeds_ceiling
        dsum_malformed
        dsum_structured_inner_unsupported
        engine_pin_backend_conflict
        fa_not_representable
        fa_rank_exceeds_dim
        fa_rank_invalid
        fa_rank_missing
        family_argument_not_recognised
        fb_cov_missing_matrix
        fb_cov_type_unknown
        fit_lacks_posterior_draws
        fixed_offset_unsupported
        fixed_smoother_not_supported
        formula_not_two_sided
        grammar_argument_duplicated
        grammar_argument_empty
        grammar_arity_too_few
        grammar_arity_too_many
        grammar_brms_known_matrices_unsupported
        grammar_brms_with_asreml_terms
        grammar_named_argument_unknown
        heterogeneous_residual_factor_not_in_cell_key
        heterogeneous_residual_family_has_no_sigma
        heterogeneous_residual_multiple_factors
        inla_gate_refused
        inla_program_failed
        inla_response_scale_unvalidated
        inla_variable_used_twice
        interaction_not_representable
        known_matrices_data_name_collision
        known_matrix_dim_mismatch
        known_matrix_dimnames_mismatch
        known_matrix_level_mismatch
        loo_requires_sampler_draws
        low_rank_rank_exceeds_basis
        low_rank_rank_invalid
        low_rank_scheme_required
        low_rank_smooth_unsupported
        met_summary_not_available
        missing_covariate_not_supported
        missing_response_refused
        model_matrix_not_recoverable
        multivariate_smooth_unsupported
        nested_index_name_collision
        nested_variable_not_in_data
        numeric_variable_in_random_interaction
        penalty_block_overlap
        penalty_block_targeting_unsupported
        penalty_candidate_columns_unscaled
        penalty_factor_contrast_ambiguous
        penalty_hyperparameter_invalid
        penalty_reference_model_inadequate
        penalty_runtime_budget_exceeded
        penalty_sampler_diagnostics_failed
        penalty_scaling_column_constant
        pol_degree_invalid
        pp_check_requires_predictive_draws
        precision_not_in_known_matrices
        precision_not_positive_definite
        precision_not_square
        precision_not_symmetric
        prior_argument_duplicated
        prior_argument_missing
        prior_argument_unknown
        prior_distribution_not_a_call
        prior_distribution_unknown
        prior_family_unsupported_inla
        prior_hyperparameter_not_scalar
        prior_hyperparameter_out_of_domain
        prior_not_translatable_for_backend
        prior_spec_empty
        prior_spec_not_formula
        prior_spec_not_two_sided
        prior_target_argument_missing
        prior_target_not_a_single_name
        prior_target_not_in_model
        prior_target_unsupported
        projection_family_unsupported
        projection_search_unsupported
        projection_structure_unsupported
        representation_unknown_for_preflight
        residual_type_unsupported_for_aggregation
        response_not_in_data
        review_code_backend_unsupported
        row_count_exceeds_integer
        selection_base_not_allowlisted
        selection_base_shape_unreadable
        selection_candidates_not_constant_within_term
        selection_covariance_not_buildable
        selection_covariance_not_positive_definite
        selection_decision_invalid
        selection_dense_reference_too_large
        selection_fallback_not_authorised
        selection_family_not_supported
        selection_fold_covariance_not_factorisable
        selection_fold_id_invalid
        selection_fold_level_absent
        selection_fold_not_usable
        selection_fold_variable_not_in_data
        selection_leakage_detected
        selection_missing_response_estimand_unnamed
        selection_missing_response_in_objective
        selection_offset_absorbed_twice
        selection_penalty_blocks_not_expressible
        selection_penalty_coordinate_basis_dependent
        selection_penalty_family_unknown
        selection_penalty_family_withdrawn
        selection_penalty_spec_malformed
        selection_penalty_target_ambiguous
        selection_penalty_target_not_in_model
        selection_penalty_target_unsupported
        selection_preprocessing_leak
        selection_provider_absent
        selection_prox_mode_not_converged
        selection_rank_deficient_protected
        selection_reduction_absorption_unsupported
        selection_reduction_level_mismatch
        selection_reduction_term_not_in_fit
        selection_resamples_underspecified
        selection_scoring_unit_mismatch
        selection_scoring_unit_undeclared
        selection_seed_required
        selection_strata_not_constant_within_unit
        selection_weights_absorbed_twice
        selection_weights_loss_route_requires_residual_only
        selection_weights_not_positive
        smooth_variable_not_in_data
        stan_cannot_represent_ar1_field
        stan_cannot_represent_ar1_residual
        stan_cannot_represent_structured_cov
        stan_cannot_represent_structured_residual
        str_not_representable
        stream_response_must_be_a_column
        tensor_smooth_unsupported
        term_in_fixed_and_random
        triangulate_incomparable_fits
        unknown_backend
        unsupported_family
        update_call_not_reconstructable
        update_unnamed_argument_not_supported
        us_gxe_prior_not_translatable_inla
        variogram_requires_design_index
        vm_redundant_specification
        vpredict_component_not_in_posterior
        vpredict_functional_not_one_value_per_draw
        vpredict_no_functional
        weights_not_aggregatable
        weights_not_supported
        weights_requires_gaussian
      
      .backend_registry: locked=TRUE n_keys=2
        brms
        inla
      

