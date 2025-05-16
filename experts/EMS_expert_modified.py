import torch
import torch.nn as nn
from typing import Optional, List
from transformers import AutoTokenizer, PreTrainedTokenizerBase
from transformers import EsmForMaskedLM
from transformers.tokenization_utils_base import BatchEncoding
from evo_prot_grad.experts.base_experts import ProteinLMExpert
import evo_prot_grad.common.embeddings as embeddings


class EsmExpertModified(ProteinLMExpert):
    """Expert wrapper for HuggingFace ESM models that outputs logits from one-hot inputs.
    Removes internal scoring (e.g., masked log-likelihood); delegates scoring to VariantScoring.
    """
    def __init__(self,
                 temperature: float,
                 scoring_strategy: str,
                 model: Optional[nn.Module] = None,
                 tokenizer: Optional[PreTrainedTokenizerBase] = None,
                 device: str = 'cpu'):
        if model is None and tokenizer is None:
            model = EsmForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D")
            tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
        elif model is None or tokenizer is None:
            raise ValueError("EsmExpert requires both `model` and `tokenizer` to be specified.")

        super().__init__(
            temperature,
            model,
            tokenizer.get_vocab(),
            scoring_strategy,
            device
        )
        self.tokenizer = tokenizer
        self.model.esm.embeddings.word_embeddings = embeddings.OneHotEmbedding(
            model.esm.embeddings.word_embeddings
        )

    def _get_last_one_hots(self) -> torch.Tensor:
        """Returns the one-hot tensors most recently passed as input."""
        return self.model.esm.embeddings.word_embeddings.one_hots

    def tokenize(self, inputs: List[str]) -> BatchEncoding:
        """Tokenize input sequences for preprocessing and sampling (no scoring)."""
        return self.tokenizer(inputs, add_special_tokens=False, return_tensors="pt").to(self.device)

    def forward(self, x_oh: torch.Tensor) -> torch.Tensor:
        """Forward pass using one-hot encoded input. Returns logits (no scoring here)."""
        return self.model(inputs_embeds=x_oh).logits
