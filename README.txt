

to use

1) need shiny, GGtools, GGdata, reactome.db, illuminaHumanv1.db,
    and all dependencies

2) start runApp() over the ui.R and server.R here

3) under 'Feature to be modeled' pick 'pathway', then "Apply Changes"

4) new controls appear for "select the pathway..."  pick the 
  "HCN channels" pathway under H-Q selector

5) optional: set Expression PCs to 10

6) Apply Changes and see the table of probes to be used

7) Set the "Configuration complete" button to "yes" and Apply Changes

8) Wait.  Eventually a new table of scores and a plot will appear.  You
can update the plot using "index for plotting" control

This does assume a two-core machine; it will not fail on one-core;
you can increase cores in use changing the RegisterDoParallel call
in server.R
