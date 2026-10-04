"""Minimal path-token Transformer for busbar centerline corrections."""

import math

import torch
from torch import nn


class PathTransformer(nn.Module):
    """Predict local normal/binormal offsets and cross-section rotation."""

    def __init__(self):
        super().__init__()
        self.input_projection = nn.Linear(13, 64)
        layer = nn.TransformerEncoderLayer(
            d_model=64, nhead=4, dim_feedforward=128, dropout=0.0,
            batch_first=True, norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(
            layer, num_layers=2, enable_nested_tensor=False,
        )
        self.head = nn.Linear(64, 3)
        # The untrained smoke run is an identity correction, not a prediction.
        nn.init.zeros_(self.head.weight)
        nn.init.zeros_(self.head.bias)

    def forward(self, features):
        if features.ndim != 3 or features.shape[-1] != 13:
            raise ValueError("features must have shape [batch, nodes, 13]")
        raw = self.head(self.encoder(self.input_projection(features)))
        offsets = 2.0 * torch.tanh(raw[..., :2])
        theta = (math.pi / 2.0) * torch.tanh(raw[..., 2:3])
        return torch.cat((offsets, theta), dim=-1)
