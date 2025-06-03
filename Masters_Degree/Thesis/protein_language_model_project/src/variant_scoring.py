import torch


class VariantScoring:
    """Every Expert has a VariantScoring object to use to score variant sequences
    containing multiple mutations.

    Supported scoring strategies:
    1) `attribute_value` - Uses a model's predicted attribute value for a given variant,
        normalized by subtracting the wildtype's predicted value.
    2) `pseudolikelihood_ratio`
    3) `mutant_marginal` 

    See: https://www.biorxiv.org/content/10.1101/2021.07.09.450648v2 for (2-3).
    """

    def __init__(self, scoring_strategy: str):
        self.scoring_strategy = scoring_strategy
        self.wt_score_cache = None

    def pseudolikelihood(self, x_oh: torch.Tensor, logits: torch.Tensor) -> torch.Tensor:
        """Compute pseudo-log-likelihood (pll) score.

        Args:
            x_oh: one-hot encoded variant sequences [parallel_chains, seq_len, vocab_size]
            logits: model logits [parallel_chains, seq_len, vocab_size]
        Returns:
            Tensor of same shape
        """
        return x_oh * torch.nn.functional.log_softmax(logits, dim=-1)

    def pseudolikelihood_ratio(self, x_oh: torch.Tensor, logits: torch.Tensor) -> torch.Tensor:
        """PLL ratio of variant vs wildtype.

        Args:
            x_oh: one-hot encoded variant sequences
            logits: predicted logits
        Returns:
            Tensor of shape [parallel_chains]
        """
        if self.wt_score_cache is None:
            raise ValueError("Wildtype pseudolikelihood must be set before calling the expert with `init_wildtype`.")

        # Debug shape prints
        print("x_oh.shape:", x_oh.shape)
        print("logits.shape:", logits.shape)
        print("pseudolikelihood(x_oh, logits).sum(dim=[1,2]).shape:",
              self.pseudolikelihood(x_oh, logits).sum(dim=[1, 2]).shape)
        print("wt_score_cache.shape:", self.wt_score_cache.shape)

        return self.pseudolikelihood(x_oh, logits).sum(dim=[1, 2]) - self.wt_score_cache.repeat(x_oh.shape[0])

    def mutant_marginal(self, x_oh: torch.Tensor, logits: torch.Tensor, wt_oh: torch.Tensor) -> torch.Tensor:
        """Mutant marginal scoring mechanism.

        Args:
            x_oh: one-hot encoded variant sequences
            logits: predicted logits
            wt_oh: one-hot encoded wild type sequence
        Returns:
            Tensor of shape [parallel_chains]
        """
        # We don't need to explicitly sum only over mutation locations,
        # because the differences at non-mutation locations are always 0
        return (self.pseudolikelihood(x_oh, logits) -
                self.pseudolikelihood(wt_oh, logits)).sum(dim=[1, 2])

    def __call__(self, x_oh: torch.Tensor, x_pred: torch.Tensor, wt_oh: torch.Tensor) -> torch.Tensor:
        """Dispatch to the appropriate scoring function.

        Args:
            x_oh: one-hot encoded variant sequence
            x_pred: model prediction (e.g., logits)
            wt_oh: one-hot encoded wildtype sequence
        Returns:
            variant_score: Tensor of shape [parallel_chains]
        """
        if self.scoring_strategy == "attribute_value":
            if self.wt_score_cache is None:
                raise ValueError("Wildtype attribute value must be set before calling the expert with `init_wildtype`.")
            return x_pred - self.wt_score_cache.repeat(x_oh.shape[0])
        elif self.scoring_strategy == "pseudolikelihood_ratio":
            return self.pseudolikelihood_ratio(x_oh, x_pred)
        elif self.scoring_strategy == "mutant_marginal":
            return self.mutant_marginal(x_oh, x_pred, wt_oh)
        else:
            raise ValueError(f"Invalid scoring strategy: {self.scoring_strategy}")

    def cache_wt_score(self, wt_oh: torch.Tensor, wt_pred: torch.Tensor) -> None:
        """Caches the score value for wildtype protein.

        Args:
            wt_oh: one-hot encoded wt protein [1, seq_len, vocab_size]
            wt_pred: model output for wt protein [1, *]
        """
        if self.scoring_strategy == "attribute_value":
            self.wt_score_cache = wt_pred
        elif self.scoring_strategy == "pseudolikelihood_ratio":
            self.wt_score_cache = self.pseudolikelihood(wt_oh, wt_pred).sum(dim=[1, 2])
