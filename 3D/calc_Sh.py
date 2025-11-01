import zipfile
import pandas as pd

timesteps_per_t = 2

with zipfile.ZipFile("data.zip") as z:
    df = pd.concat([pd.read_csv(z.open(f"data.{i}.csv")) for i in range(20*timesteps_per_t, 30*timesteps_per_t)],
                    ignore_index=True, sort=False)
    print(f"Sh (t=20-30 average) = {df['Sh'].mean():.1f}")

    df = pd.concat([pd.read_csv(z.open(f"data.{i}.csv")) for i in range(50*timesteps_per_t, 100*timesteps_per_t)],
                    ignore_index=True, sort=False)
    print(f"Sh (t=50-100 average) = {df['Sh'].mean():.1f}")