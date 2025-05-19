import torch
import torch.nn as nn
from typing import Optional, List
from transformers import AutoTokenizer, PreTrainedTokenizerBase
from transformers import EsmForMaskedLM
from transformers.tokenization_utils_base import BatchEncoding
from evo_prot_grad.experts.base_experts import ProteinLMExpert
import evo_prot_grad.common.embeddings as embeddings


class EsmExpertModified(ProteinLMExpert):
    """Expert baseclass for HuggingFace protein language models from the ESM family.
    Implements abstract methods `_get_last_one_hots` and `tokenize`.
    Swaps out the `EsmForMaskedLM.esm.embeddings.word_embeddings` layer
    for a `evo_prot_grad.common.embeddings.OneHotEmbedding` layer.
    In the future, will have a masking objective added as well, to prevent biases towards certain amino acids
    
    """
    def __init__(self,
                 temperature: float,
                 scoring_strategy: str,
                 model: Optional[nn.Module] = None,
                 tokenizer: Optional[PreTrainedTokenizerBase] = None,
                 device: str = 'cpu'):
        if model is None and tokenizer is None:
            model = EsmForMaskedLM.from_pretrained("facebook/esm2_t33_650M_UR50D")
            tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")
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
