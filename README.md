# tracedsp - miniSEED time series processesor

A general purpose tool to apply signal processing to miniSEED data.

## Documentation

For usage infromation see the [tracedsp manual](doc/tracedsp.md) in the
'doc' directory.

## Downloading and building

The [releases](https://github.com/earthscope/tracedsp/releases) area
contains release versions.

In most environments a simple 'make' will build the program.

The CC and CFLAGS environment variables can be used to configure
the build parameters.

## Caveats

This program has lots of signal processing, selection, etc. capability
but not a lot of guardrails for illogical combinations of options.
The caller is encouraged to test specific combinations before applying
to research data.

As to development, the functionality has outgrown the original
code design.  It works, but it's not beautiful code.

## License

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Copyright (C) 2023 Chad Trabant, EarthScope Data Service
