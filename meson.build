project('bore', 'c', 'cpp',
  version : '1',
  default_options : ['c_std=c17', 'cpp_std=c++17', 'warning_level=3', 'buildtype=release', 'b_ndebug=if-release'])

static_link_gsl = get_option('static_gsl')

sources = [
     'src/bore.c'
]

vapoursynth_dep = dependency('vapoursynth', version: '>=55').partial_dependency(compile_args : true, includes : true)

gsl_dep = dependency('gsl', static : static_link_gsl)

deps = [vapoursynth_dep, gsl_dep]

shared_module('bore', sources,
  dependencies : deps,
  install : true,
  install_dir : join_paths(vapoursynth_dep.get_pkgconfig_variable('libdir'), 'vapoursynth'),
  gnu_symbol_visibility : 'hidden'
)
