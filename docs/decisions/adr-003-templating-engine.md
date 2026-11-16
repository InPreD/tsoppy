# Templating engine

## Context and Problem Statement

The aim here is to separate text generating functionality from the rest of the code so that:
1) all the text outputs are kept in a well defined and easy to find place,
1) it is easy to edit the text when needed,
2) data structures generated with python code and containing output information can be easily rendered.

Select a templating engine to dynamically generate output texts, consider readability of the template code for others as well as how easy it is to modify the code.
Considering [jinja2](https://jinja.palletsprojects.com/en/stable/) and [mako](https://www.makotemplates.org/).

Both engines were used to create a template for the report_predispositions module to see how they work in the real world example. In this case, the templates had to render some string variables, floats, join content of a string list, work with a dictionary and a data frame containing simple string fields (a more complex data frame needs to be transformed into such a simple one to explicitly state what should be rendered and how). This probably covers most data types (if not all) that we want to render into the output files.
<!--Include all relevant information necessary for making a decision here.-->

## Decision

From the two tested possibilities, `mako` was much better in enforcing proper formatting and rendering of various data types (int, float, complex data frames, etc.). I see this as a positive trait. `jinja` did not complain much, rendered something and usually it looked ok, but debugging or changing default behavior seem to be taking more time in the end.
Thus `mako` seems to be a better option going forward and we will start using that package.

<!--State what was decided.-->

## Consequences

- create template folder in each module where templates will be used
- set the absolute path for the template folder inside the container (see the code example below)
- install mako package in the repo
- use the package to work with the templates

<!--- List general consequences.-->

### Positive

- separation of the text templates to a dedicated place 
- ease to modify and read the text 
<!--- List positive consequences.-->

### Negative

- learning curve
<!--- List negative consequences.-->

### Neutral

<!--- List neutral consequences.-->

## Code example 

In the python script:
```
from mako.template import Template as makoTemplate

my_template = makoTemplate(
        filename='/workspaces/tsoppy/src/tsoppy/<MODULE>/templates/<TEMPLATE_NAME>.mako')
        
output = my_template.render(str_variable = val1, float_variable = val2, list_variable = val3, dataframe_variable = val4)
```

In the template:
```
This is an example of a mako template. 
Here is an input string value: ${str_variable}
When formatting floats, trailing zeros can be kept using {:.4f}.format(variable): ${ "{:.4f}".format(float_variable) }  
One can print all the list elements using join function: ${", ".join(list_variable)}
And to work with a simple, string containing data frame, one can do:
Print header:
${ "\t".join(dataframe_variable.columns) }
Iterate through all the rows:
% for row in dataframe_variable.iter_rows():
${ "\t".join(row) }
% endfor
```

For a bit more, have a look [here](https://github.com/InPreD/tsoppy/blob/report_predispositions/src/tsoppy/report_predispositions/templates/predispositions.mako)
or dive into the [documentation](https://www.makotemplates.org/)
