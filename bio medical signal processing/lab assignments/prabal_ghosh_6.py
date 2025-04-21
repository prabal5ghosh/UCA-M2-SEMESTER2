import numpy as np
import scipy.io as sio
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset, random_split
from sklearn.metrics import mean_squared_error

# === LOAD & PREP DATA ===

# Load Ra files
ra1 = sio.loadmat('C:\\Users\\praba\\Documents\\GitHub\\UCA-M2-SEMESTER2\\bio medical signal processing\\lab assignments\\Atrial_prabal_ghosh\\Ra1.mat')['Ra1']
ra2 = sio.loadmat('C:\\Users\\praba\\Documents\\GitHub\\UCA-M2-SEMESTER2\\bio medical signal processing\\lab assignments\\Atrial_prabal_ghosh\\Ra2.mat')['Ra2']
ra3 = sio.loadmat('C:\\Users\\praba\\Documents\\GitHub\\UCA-M2-SEMESTER2\\bio medical signal processing\\lab assignments\\Atrial_prabal_ghosh\\Ra3.mat')['Ra3']

# Merge into one big 3D array → shape: (12, 15000, 75)
Xa_all = np.concatenate([ra1, ra2, ra3], axis=2)

# Extract only V1–V6 (indices 6 to 11)
Xa_v1_v6 = Xa_all[6:12, :, :]  # shape: (6, 15000, 75)

# Rearrange to (75 patients, 15000 time steps, 6 leads)
X = np.transpose(Xa_v1_v6, (2, 1, 0))  # shape: (75, 15000, 6)

# Load Rva files
rva1 = sio.loadmat('C:\\Users\\praba\\Documents\\GitHub\\UCA-M2-SEMESTER2\\bio medical signal processing\\lab assignments\\Atrial_prabal_ghosh\\Rva1.mat')['Rva1']
rva2 = sio.loadmat('C:\\Users\\praba\\Documents\\GitHub\\UCA-M2-SEMESTER2\\bio medical signal processing\\lab assignments\\Atrial_prabal_ghosh\\Rva2.mat')['Rva2']
rva3 = sio.loadmat('C:\\Users\\praba\\Documents\\GitHub\\UCA-M2-SEMESTER2\\bio medical signal processing\\lab assignments\\Atrial_prabal_ghosh\\Rva3.mat')['Rva3']

Rva_all = np.concatenate([rva1, rva2, rva3], axis=2)
Rva_v1_v6 = Rva_all[6:12, :, :]  # Extract V1–V6
Rva = np.transpose(Rva_v1_v6, (2, 1, 0))  # Rearrange to (75, 15000, 6)

# Convert data to PyTorch tensors
X_tensor = torch.tensor(X, dtype=torch.float32)
Rva_tensor = torch.tensor(Rva, dtype=torch.float32)

# Create a dataset and split into train and test sets
dataset = TensorDataset(Rva_tensor, X_tensor)
train_size = int(0.8 * len(dataset))
test_size = len(dataset) - train_size
train_dataset, test_dataset = random_split(dataset, [train_size, test_size])

train_loader = DataLoader(train_dataset, batch_size=4, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=4, shuffle=False)

# === DEFINE BiLSTM MODEL ===
class BiLSTMModel(nn.Module):
    def __init__(self):
        super(BiLSTMModel, self).__init__()
        self.lstm1 = nn.LSTM(input_size=6, hidden_size=64, num_layers=1, batch_first=True, bidirectional=True)
        self.dropout1 = nn.Dropout(0.2)
        self.lstm2 = nn.LSTM(input_size=128, hidden_size=64, num_layers=1, batch_first=True, bidirectional=True)
        self.dropout2 = nn.Dropout(0.2)
        self.fc = nn.Linear(128, 6)

    def forward(self, x):
        x, _ = self.lstm1(x)
        x = self.dropout1(x)
        x, _ = self.lstm2(x)
        x = self.dropout2(x)
        x = self.fc(x)
        return x

# Initialize the model, loss function, and optimizer
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = BiLSTMModel().to(device)
criterion = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

# === TRAIN THE MODEL ===
epochs = 10
for epoch in range(epochs):
    model.train()
    train_loss = 0.0
    for inputs, targets in train_loader:
        inputs, targets = inputs.to(device), targets.to(device)

        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, targets)
        loss.backward()
        optimizer.step()

        train_loss += loss.item()

    print(f"Epoch {epoch + 1}/{epochs}, Loss: {train_loss / len(train_loader):.4f}")

# === EVALUATE THE MODEL ===
model.eval()
y_true = []
y_pred = []
with torch.no_grad():
    for inputs, targets in test_loader:
        inputs, targets = inputs.to(device), targets.to(device)
        outputs = model(inputs)
        y_true.append(targets.cpu().numpy())
        y_pred.append(outputs.cpu().numpy())

y_true = np.concatenate(y_true, axis=0)
y_pred = np.concatenate(y_pred, axis=0)

# Calculate Mean Squared Error
mse = mean_squared_error(y_true.flatten(), y_pred.flatten())
print("Test MSE:", mse)