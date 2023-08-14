#!/usr/bin/env python3

# Modules for vanilla variational autoencoder

import torch
import torch.nn.functional as F
from torch import nn, optim
from utils import EarlyStopping

class VAE(torch.nn.Module):
    
    def __init__(self, mask_tf, depth, **kwargs):
        """ Constructor for class pathway VAE """
        super(VAE, self).__init__()
        self.mask_tf = torch.BoolTensor(mask_tf)
        self.n_genes = len(self.mask_tf)
        self.n_tfs = int(self.mask_tf.long().sum())
        self.depth = depth
        self.init_w = kwargs.get('init_w', True)
        self.beta = kwargs.get('beta', 0.05)
        self.alpha = kwargs.get('alpha', 0.5)
        self.dropout = kwargs.get('dropout', 0.3)
        self.path_model = kwargs.get('path_model', "trained_vae.pt")
        self.dev = kwargs.get('device', torch.device('cpu'))
        
        encoder = []
        for i in range(self.depth):
            encoder.append(nn.Linear(self.n_genes, self.n_genes))
            encoder.append(nn.BatchNorm1d(self.n_genes))
            encoder.append(nn.ReLU())
            encoder.append(nn.Dropout(self.dropout))   
        self.encoder = nn.Sequential(*encoder)
        
        self.mean = nn.Sequential(nn.Linear(self.n_genes, self.n_tfs),
                                    nn.Dropout(self.dropout))
        self.logvar = nn.Sequential(nn.Linear(self.n_genes, self.n_tfs),
                                    nn.Dropout(self.dropout))

        decoder = []
        decoder.append(nn.Linear(self.n_tfs, self.n_genes))
        for i in range(self.depth - 1):
            decoder.append(nn.BatchNorm1d(self.n_genes))
            decoder.append(nn.ReLU())
            decoder.append(nn.Dropout(self.dropout))
            decoder.append(nn.Linear(self.n_genes, self.n_genes))
        decoder.append(nn.ReLU())#force to be non-negative
        self.decoder = nn.Sequential(*decoder)
        
        if self.init_w:
            self.encoder.apply(self._init_weights)
            self.decoder.apply(self._init_weights)           
        
    def _init_weights(self, layer):
        """ Initialize weights of layer with Xavier uniform"""
        if type(layer)==nn.Linear:
            torch.nn.init.xavier_uniform_(layer.weight)
        return
    
    def sample_latent(self, mu, logvar):
        """ Sample latent space with reparametrization trick. First convert to std, sample normal(0,1) and get Z."""
        std = logvar.mul(0.5).exp_()
        eps = torch.FloatTensor(std.size()).normal_().to(self.dev)
        eps = eps.mul_(std).add_(mu)
        return eps
    
    def encode(self, X):
        """ Encode data """
        y = self.encoder(X)
        mu, logvar = self.mean(y), self.logvar(y)
        z = self.sample_latent(mu, logvar)
        z = F.relu(z)#force to be non-negative
        return z, mu, logvar
    
    def decode(self, z):
        """ Decode data """
        X_rec = self.decoder(z)
        return X_rec
    
    def forward(self, X):
        """ Forward pass through full network"""
        z, mu, logvar = self.encode(X)
        X_rec = self.decode(z)
        return X_rec, z, mu, logvar

    def vae_loss(self, X, X_rec, Z, Z_rec, mu, logvar):
        """ Custom loss for VAE """
        kld = -0.5 * torch.sum(1. + logvar - mu.pow(2) - logvar.exp(), )
        loss_1 = (1 - self.beta) * F.mse_loss(X_rec, X, reduction="sum") + self.beta * kld
        loss_2 = F.mse_loss(Z_rec, Z, reduction="sum")
        return torch.mean((1 - self.alpha) * loss_1 + self.alpha * loss_2)

    def train_model(self, train_loader, learning_rate, n_epochs, train_patience, test_patience, test_loader=False, save_model=True):
        """ Train VAE """
        epoch_hist = {}
        epoch_hist['train_loss'] = []
        epoch_hist['valid_loss'] = []
        optimizer = optim.Adam(self.parameters(), lr=learning_rate, weight_decay=5e-4)
        train_ES = EarlyStopping(patience=train_patience, verbose=True, mode='train')
        if test_loader:
            valid_ES = EarlyStopping(patience=test_patience, verbose=True, mode='valid')
        
        # Train
        for epoch in range(n_epochs):
            loss_value = 0
            self.train()
            for train in train_loader:
                train = train.to(self.dev)
                x_train = train
                z_train = train[:, self.mask_tf]
                optimizer.zero_grad()
                x_rec, z_rec, mu, logvar = self.forward(x_train)
                loss = self.vae_loss(x_train, x_rec, z_train, z_rec, mu, logvar)
                loss_value += loss.item()
                loss.backward()
                optimizer.step()

            # Get epoch loss
            epoch_loss = loss_value / (len(train_loader) * train_loader.batch_size)
            epoch_hist['train_loss'].append(epoch_loss)
            train_ES(epoch_loss)
            # Eval
            if test_loader:
                self.eval()
                test_dict = self.test_model(test_loader)
                test_loss = test_dict['loss']
                epoch_hist['valid_loss'].append(test_loss)
                valid_ES(test_loss)
                print('[Epoch %d] | loss: %.3f | test_loss: %.3f |'%(epoch+1, epoch_loss, test_loss), flush=True)
                if valid_ES.early_stop or train_ES.early_stop:
                    print('[Epoch %d] Early stopping' % (epoch+1), flush=True)
                    break
            else:
                print('[Epoch %d] | loss: %.3f |' % (epoch + 1, epoch_loss), flush=True)
                if train_ES.early_stop:
                    print('[Epoch %d] Early stopping' % (epoch+1), flush=True)
                    break
       
        # Save model
        if save_model:
            print('Saving model to ...', self.path_model)
            torch.save(self.state_dict(), self.path_model)

        return epoch_hist

    def test_model(self, loader):
        """Test model on input loader."""
        test_dict = {}
        loss = 0
        loss_func = self.vae_loss
        self.eval()
        with torch.no_grad():
            for data in loader:
                data = data.to(self.dev)
                x = data
                z = data[:, self.mask_tf]
                x_rec, z_rec, mu, logvar = self.forward(x)
                loss += loss_func(x, x_rec, z, z_rec, mu, logvar).item()
        test_dict['loss'] = loss/(len(loader)*loader.batch_size)
        return test_dict




