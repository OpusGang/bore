project('bore', 'c', 'cpp',
  version : '1.0.0',
  default_options : ['c_std=c17', 'cpp_std=c++20', 'warning_level=3', 'buildtype=release', 'b_ndebug=if-release'])

static_link_gsl = get_option('static_gsl')

python = find_program('python3', 'python', required: true)

includedirs = ['src']

sources = [
     'src/common.c',
     'src/common.h'
]

deps = []

libtype = get_option('libtype')

vs_dep_for_paths = dependency('vapoursynth', version: '>=55', required: false)
avs_dep_for_paths = dependency('avisynth', required: false)

if libtype in ['vapoursynth', 'both']
    sources += ['src/bore.c']

    vs_header_dep = vs_dep_for_paths.partial_dependency(compile_args: true, includes: true)

    if vs_header_dep.found()
        message('VapourSynth headers found via dependency system. Using its include paths.')
        deps += [vs_header_dep]
    else
        message('VapourSynth headers NOT found via dependency system.')
        message('If building VapourSynth support, ensure VapourSynth headers are in the include path (e.g., via CFLAGS or INCLUDE env var).')
    endif
endif

if libtype in ['avisynth', 'both']
    sources += [
        'avs_c_api_loader/avs_c_api_loader.cpp',
        'avs_c_api_loader/avs_c_api_loader.hpp',
        'src/avs_bore.cpp'
    ]

    avs_header_dep = avs_dep_for_paths.partial_dependency(compile_args: true, includes: true)

    if avs_header_dep.found()
        message('AviSynthPlus headers found via dependency system. Using its include paths.')
        deps += [avs_header_dep]
    else
        message('AviSynthPlus headers NOT found via dependency system.')
        message('If building AviSynthPlus support, ensure AviSynthPlus headers are in the include path (e.g., via CXXFLAGS or INCLUDE env var).')
    endif
endif

# Windows Version Information
if host_machine.system() == 'windows'
    windows = import('windows')

    full_project_version_str = meson.project_version() # e.g., '1.0.0-alpha', '1.0', '1'

    version_parts_split_by_hyphen = full_project_version_str.split('-')
    core_version_str = version_parts_split_by_hyphen[0] # e.g., '1.0.0', '1.0', '1'

    project_version_parts = core_version_str.split('.')
    version_major = project_version_parts[0].to_int()
    version_minor = project_version_parts.length() > 1 ? project_version_parts[1].to_int() : 0
    version_patch = project_version_parts.length() > 2 ? project_version_parts[2].to_int() : 0
    version_build = project_version_parts.length() > 3 ? project_version_parts[3].to_int() : 0

    conf_data = configuration_data()
    conf_data.set('VERSION_MAJOR', version_major)
    conf_data.set('VERSION_MINOR', version_minor)
    conf_data.set('VERSION_PATCH', version_patch)
    conf_data.set('VERSION_BUILD', version_build)
    conf_data.set_quoted('PROJECT_VERSION_STRING_FULL', full_project_version_str)

    target_apis_str = ''
    if libtype == 'vapoursynth'
        target_apis_str = 'VapourSynth'
    elif libtype == 'avisynth'
        target_apis_str = 'AviSynthPlus'
    else
        target_apis_str = 'VapourSynth and AviSynthPlus'
    endif

    conf_data.set_quoted('FILE_DESCRIPTION', 'Border fixing filter for ' + target_apis_str)
    conf_data.set_quoted('INTERNAL_NAME', meson.project_name().to_lower())
    conf_data.set_quoted('ORIGINAL_FILENAME', 'lib' + meson.project_name().to_lower() + '.dll')
    conf_data.set_quoted('PRODUCT_NAME', meson.project_name())

    is_debug_build = get_option('buildtype') == 'debug'
    is_prerelease = version_parts_split_by_hyphen.length() > 1 ? true : false

    conf_data.set('PROJECT_DEBUG_BUILD', is_debug_build ? 1 : 0)
    conf_data.set('PROJECT_IS_PRERELEASE', is_prerelease ? 1 : 0)

    windows_rc_configured = configure_file(
        input: 'src/bore.rc.in',
        output: 'bore.rc',
        configuration: conf_data
    )

    compiled_rc = windows.compile_resources(
        windows_rc_configured,
        include_directories: includedirs
    )
    sources += compiled_rc
endif

gsl_dep = dependency('gsl', static : static_link_gsl)

deps += gsl_dep

bore_module = shared_module('bore', sources,
  dependencies : deps,
  install : false,
  link_args : '-s',
  gnu_symbol_visibility : 'hidden'
)

# Manual installation
install_dirs_to_use = []

if libtype in ['vapoursynth', 'both']
    vs_installdir = ''
    if host_machine.system() == 'windows'
        vs_installdir = join_paths(get_option('libdir'), 'vapoursynth')
    else
        vs_pkg_libdir = vs_dep_for_paths.get_variable(pkgconfig: 'libdir', default_value: '')
        if vs_pkg_libdir == ''
            message('Warning: VapourSynth pkg-config "libdir" not found. Using default libdir for VapourSynth plugin installation.')
            vs_installdir = get_option('prefix') / get_option('libdir') / 'vapoursynth'
        else
            vs_installdir = join_paths(vs_pkg_libdir, 'vapoursynth')
        endif
    endif
    install_dirs_to_use += [vs_installdir]
endif

if libtype in ['avisynth', 'both']
    avs_installdir = ''
    if host_machine.system() == 'windows'
        avs_installdir = join_paths(get_option('libdir'), 'avisynth')
    else
        avs_pkg_libdir = avs_dep_for_paths.get_variable(pkgconfig: 'libdir', default_value: '')
        if avs_pkg_libdir == ''
            message('Warning: AviSynthPlus pkg-config "libdir" not found. Using default libdir for AviSynthPlus plugin installation.')
            avs_installdir = get_option('prefix') / get_option('libdir') / 'avisynth'
        else
            avs_installdir = join_paths(avs_pkg_libdir, 'avisynth')
            message('Using Avisynth pkg-config libdir for plugin installation: ' + avs_installdir)
        endif
    endif
    install_dirs_to_use += [avs_installdir]
endif

if install_dirs_to_use.length() == 0 and libtype != 'none'
    error('No installation directory determined for libtype: ' + libtype)
endif

if install_dirs_to_use.length() > 0
    message('Setting up installation of bore module to: ' + ' | '.join(install_dirs_to_use))

    install_script_path = join_paths(meson.current_source_dir(), 'install_plugin.py')

    foreach installdir : install_dirs_to_use
        meson.add_install_script(
            python.full_path(),
            install_script_path,
            bore_module,
            installdir
        )
    endforeach
else
    if libtype != 'none'
        message('bore module built but not installed as no valid install directories were determined for libtype: ' + libtype)
    endif
endif
