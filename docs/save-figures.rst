Save figures
============

This is a short guide on how to save high-resolution figures from the dynophore notebook.

NGLView figures
---------------

You can download the ``nglview`` view as an image as follows:

.. code-block:: python

    view = dyno.view3d.show(..)
    view


Adjust the molecule as you wish it to be saved. Then save it!

.. code-block:: python

    view.download_image(filename='screenshot.png', factor=4, antialias=True, trim=False, transparent=False)


Find more details in the ``nglview`` 
`documentation <http://nglviewer.org/nglview/latest/_api/nglview.widget.html#nglview.widget.NGLWidget.download_image>`_.

Plots
-----

In the dynophore notebook, all plots are interactive thanks to 
the ``dynophores.plot.interactive`` module.

If you want to save figures, you can use the ``ipympl`` widget (disc icon); 
find more details `here <https://github.com/matplotlib/ipympl>`_. 
You can also enlarge the figure first and save it then to increase resolution.

If you want to manipulate figures or save them at a defined high resolution, 
you can use the ``dynophores.plot.static`` module, 
which follows the same signature as its interactive counterpart. 
Functions from this module will return the figure and axes handles, allowing to apply changes.

The interactive mode:

.. code-block:: python

    dyno.plot.interactive.envpartners_all_in_one(dynophore)
    # Save figure via the ipympl widget menu

The static mode:

.. code-block:: python

    fig, ax = dyno.plot.static.envpartners_all_in_one(dynophore)
    fig.savefig("myfigure.png", dpi=300)

Note that in the static mode you need to select the exact data you would like to display 
(in the interactive mode this is by interacting with the widgets).
You can check your options by typing ``dyno.plot.static.envpartners_all_in_one?`` 
followed by a tab and enter.