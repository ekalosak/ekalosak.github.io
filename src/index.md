# Everything's a graph

Welcome to my blog.

# Use panda3d instead of pygame in openai/gym
Wed 28, 2022

```python
class MyBase(panda3d.ShowBase):
  ...

class Env:
  def render():
    if mode == human:
      self.window = My
```
