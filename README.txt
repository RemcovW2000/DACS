This is the composite wing analysis toolbox I'm making. It's currently definitely a work in progress but definitely a
project I'm proud of.

The easiest way to run an example of what the code actually does, is to run the streamlit app with the following command:

streamlit run Toolbox\streamlit_app.py

In the streamlit file under /Toolbox/streamlit_app.py, the wing is fully defined, the definition is unfortunately quite
complex and difficult to operate. More complete documentation would be necessary for it to be easily usable for other
users.

I wrote some documentation which shows the overall flow of the code, as well as the class structure. It's a very simplified
overview. It's a html page made using sphinx, please check it out in /docs/build/html/index.html. In pycharm it's possible
to open the page directly from the editor. Otherwise its possible to launch from the terminal. Navigate to /docs/build/html
and execute 'index.html' as a command.

Required packages:
Numpy
streamlit
pandas
matplotlib
copy
sys
os
scipy
tqdm
math
