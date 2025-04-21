import torch
import torch.nn as nn
import torch.nn.functional as F
from layers import GraphAttentionLayer, SpGraphAttentionLayer


class multiGATModelAE(nn.Module):
    def __init__(self, nfeat, nhid, nclass, dropout, alpha, nheads, npatient):
        super(multiGATModelAE, self).__init__()
        self.dropout = dropout
                
        self.RNASeq = GraphGAT(nfeat, nhid, nclass, dropout, alpha, nheads, npatient)
        
        self.DNAm = GraphGAT(nfeat, nhid, nclass, dropout, alpha, nheads, npatient)
        
        self.CNA = GraphGAT(nfeat, nhid, nclass, dropout, alpha, nheads, npatient)

        self.dc = InnerProductDecoder(dropout, act = lambda x: x)
        
    def forward(self, rnaseq, dnam, cna, adj):
        x_r = self.RNASeq(rnaseq, adj)
        x_d = self.DNAm(dnam, adj)
        x_c = self.CNA(cna, adj)
        
        z = (x_r+ x_d+ x_c) / 3
        
        return self.dc(z), z

class InnerProductDecoder(nn.Module):
    """Decoder for using inner product for prediction."""

    def __init__(self, dropout, act=torch.sigmoid):
        super(InnerProductDecoder, self).__init__()
        self.dropout = dropout
        self.act = act

    def forward(self, z):
        z = F.dropout(z, self.dropout, training=self.training)
        adj = self.act(torch.mm(z, z.t()))
        return adj


class GraphGAT(nn.Module):
    def __init__(self, nfeat, nhid, nclass, dropout, alpha, nheads, npatient):
        super(GraphGAT, self).__init__()
        self.dropout = dropout

        self.attentions = [GraphAttentionLayer(nfeat, nhid, dropout=dropout, alpha=alpha, concat=True) for _ in range(nheads)]
        for i, attention in enumerate(self.attentions):
            self.add_module('attention_{}'.format(i), attention)

        self.out_att = GraphAttentionLayer(nhid * nheads, nclass, dropout=dropout, alpha=alpha, concat=False)
        
        attention_size = 16
        self.Wz = nn.Parameter(torch.empty(size=(nclass, attention_size)))
        nn.init.xavier_uniform_(self.Wz.data, gain=1.414) #xavier初始化
        self.Wa = nn.Parameter(torch.empty(size=(npatient, attention_size)))
        nn.init.xavier_uniform_(self.Wa.data, gain=1.414) #xavier初始化
        self.v = nn.Parameter(torch.empty(size=(attention_size, 1)))
        nn.init.xavier_uniform_(self.v.data, gain=1.414)

    def forward(self, x, adj):
        z = F.dropout(x, self.dropout, training=self.training)
        z = torch.cat([att(z, adj) for att in self.attentions], dim=1)
        z = F.dropout(z, self.dropout, training=self.training)
        z = F.elu(self.out_att(z, adj))
        a = F.tanh(torch.mm(z, Wz) + torch.mm(adj, Wa))
        a = torch.mm(a, self.v)
        a = F.softmax(a, dim = 1)
        a = a + a.T
        z = torch.mm(a, z)              
        return z





