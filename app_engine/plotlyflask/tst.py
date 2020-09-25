import seaborn as sns
# sns.set_theme(style="whitegrid")
tips = sns.load_dataset("tips")
ax = sns.swarmplot(x=tips["total_bill"])

import plotly.express as px
import plotly.io as pio

pio.renderers.default="chrome"


df = px.data.tips()
fig = px.violin(df, y="tip", x="smoker", color="sex", box=True, points="all",
          hover_data=df.columns)
fig.show()

fig = px.strip(df, x='day', y='tip')
fig.show()
