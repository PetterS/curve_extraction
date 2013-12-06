% This script compiles files only if needed
%
% It searches for the Spii library and the Vessel library in the default
% locations.
%
% On Windows, it will need help to find Eigen includes. Use the same
% version of Eigen as when compiling Spii.
%
function compile(base_path, base_name, sources, extra_args)

if nargin < 3
	sources = {};
end

if nargin < 4
	extra_args = {};
end

% Windows note: MatLabs recommended compiler:
% (Microsoft Software Development Kit (SDK) 7.1)
% Does not have support for OpenMP, see
% http://openmp.org/wp/openmp-compilers/.
%
% For Matlab 2013a:
% Visual Studio Pro 2013  has OpenMP support:
% Note visual studio version support depends on your Matlab version. 
%
% --
use_openmp = true;
if ~ispc
    CXXFLAGS = '-std=c++0x';
    
    if (use_openmp)
        extra_args{end+1} = '-lgomp';
        CXXFLAGS = [CXXFLAGS '-fopenmp'];
    end

    extra_args{end+1} = ['CXXFLAGS="\$CXXFLAGS ' CXXFLAGS '"'];

elseif (use_openmp)
    extra_args{end+1} = 'COMPFLAGS="$COMPFLAGS /openmp"'; 
end

if (use_openmp)
    extra_args{end+1} = '-DUSE_OPENMP=1';
end

% Compile script location
current_dir = fileparts(mfilename('fullpath'));

if ispc
    % Windows settings.
    
    %%%%%%%%%
    % ACTION REQUIRED:
    % Point this to your Eigen include directory.
    extra_args{end+1} = ['-I' cd '/../../../../../Programming/spii/thirdparty/Eigen/'];
    %%%%%%%%

    % Default installation directories for Spii and curve_extraction libraries.
    extra_args{end+1} = '-I"C:\Program Files\curve_extraction\include"';
    extra_args{end+1} = '-I"C:\Program Files\SPII\include"';
    extra_args{end+1} = '-L"C:\Program Files\curve_extraction\lib"';
    extra_args{end+1} = '-L"C:\Program Files\SPII\lib"';
else
    % Linux settings.
    extra_args{end+1} = '-I/usr/local/include/spii-thirdparty';
    extra_args{end+1} = '-I/usr/local/include/eigen3';
    extra_args{end+1} = '-I/usr/local/include/spii';
    extra_args{end+1} = '-I/usr/local/include/curve_extraction';
end

% Link with external libaries spii and curve_extraction.
extra_args{end+1} = '-lspii';
extra_args{end+1} = '-lcurve_extraction';


m_file_name   = [base_path filesep base_name '.m'];
mex_file_name = [base_path filesep base_name '_mex.' mexext];
cpp_file_name = [base_path filesep base_name '_mex.cpp'];

% Cpp file not found returning
if (isempty(cpp_file_name))
    fprintf('Warning: C++ file not found \n');
    return;
end


m_file   = dir(m_file_name);
mex_file = dir(mex_file_name);
cpp_file = dir(cpp_file_name);

if length(mex_file) == 1
    mex_modified = mex_file.datenum;
else
    mex_modified = 0;
end

m_modified   = m_file.datenum;
cpp_modified = cpp_file.datenum;

% If modified or not existant compile
compile_file = false;
if ~exist(mex_file_name, 'file')
    disp(mex_file_name)
        disp('Mex file not found; compiling...');
        compile_file = true;
 elseif mex_modified < cpp_modified
         disp('C++ file modfied later than Mex file; recompiling...');
         compile_file = true;
 elseif mex_modified < m_modified
         disp('M-file modfied later than Mex file; recompiling...');
         compile_file = true;
end

%% Checking additional c++ files
include_folders = {};
for i = 1 : length(sources);
	include_folders{i} = ['-I' base_path filesep fileparts(sources{i}) filesep];
	sources{i} = [base_path filesep sources{i}];
end

% Check all dependend files
for i = 1 : length(sources)
	cpp_file = dir(sources{i});
	cpp_modified = cpp_file.datenum;
	if mex_modified < cpp_modified
		compile_file = true;
	end
end

%% Compile
if compile_file
    disp(['Compiling ' cpp_file_name '...']);

    mex_argument = {'mex', ...
        cpp_file_name,  ...
        '-outdir', ...
        base_path, ...
        '-largeArrayDims', ...
        extra_args{:},...
        include_folders{:}, ...
        sources{:}};
    
    % Spaceing needed
    for i = 1:numel(mex_argument)
        mex_argument{i} = [mex_argument{i} ' '];
    end
    
    % Using eval to allow for arguments changeing c and c++ flags
    eval([mex_argument{:}])

    disp('done.');
end

