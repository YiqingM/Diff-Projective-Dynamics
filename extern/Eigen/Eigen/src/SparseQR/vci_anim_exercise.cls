%!TEX root = ./main.tex
% exam.cls
% Copyright (c) 2016 Dan A. Koschier
%
% A LaTeX2e document class for preparing exams.

\def\fileversion{0.1}
\def\filedate{2016/08/22}

\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{vci_anim_exercise}[2015/11/21 VCI Animation Exercise Template]%
%\LoadClass{scrartcl}%
\LoadClass{article}%

%\RequirePackage{titlesec}%
\RequirePackage{amsmath}
\RequirePackage{amsfonts}
\RequirePackage{amssymb}
\RequirePackage[includehead,includefoot]{geometry}
\RequirePackage[explicit, compact]{titlesec}
\RequirePackage{titling}
\RequirePackage[hidelinks]{hyperref}
\usepackage{watermark}

\RequirePackage[singlespacing]{setspace}
\RequirePackage[defaultsans]{droidsans}
\RequirePackage{enumitem}

\usepackage{ifthen}

\def\shortyear#1{\expandafter\shortyearhelper#1}
\def\shortyearhelper#1#2#3#4{#3#4}

\newboolean{@drafttrue}
\setboolean{@drafttrue}{false}
\DeclareOption{draft}{\setboolean{@drafttrue}{true}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\ProcessOptions

\setlength{\hoffset}{-0.55in}
\setlength{\voffset}{-1.0in}
\setlength{\textwidth}{500pt}
\setlength{\textheight}{680pt}
%\usepackage{showframe}


%\setlength{\voffset}{-0.75in}
%\setlength{\textheight}{650pt}
%\setlength{\marginparsep}{0pt}
%\setlength{\marginparwidth}{0pt}
%\setlength{\textwidth}{450pt}
%\setlength{\hoffset}{0pt}

%\setlength{\evensidemargin}{0.122803in}
%\setlength{\oddsidemargin}{-0.081343in}
\let\endtitlepage\relax

\titleformat{name=\section}[block]
{\bfseries\large}
{Assignment \thesection\hspace{7pt}-}
{7pt}
{#1}

\titleformat{name=\section,numberless}[block]
{\bfseries\large}
{}
{0pt}
{#1}

\titleformat{name=\subsection,numberless}[block]
{\bfseries\large}
{}
{0pt}
{\noindent#1}

\renewcommand*\familydefault{\sfdefault}
\usepackage[T1]{fontenc}

\setlength{\droptitle}{-8.2em}   % This is your set screw
\pretitle{\begin{flushleft}\bfseries\LARGE}
	\posttitle{\par\end{flushleft}}
\preauthor{\begin{flushleft}\large}
	\postauthor{\end{flushleft}}
\predate{\begin{flushleft}}
	\postdate{\par\end{flushleft}}

\titlespacing*{\section}
{0pt}{5.5ex plus 1ex minus .2ex}{3.0ex plus .2ex}

\thiswatermark{\put(235,-60){\includegraphics[scale=1.0]{logo.pdf}}}

\setlist[enumerate,1]{label={\alph*)}}
