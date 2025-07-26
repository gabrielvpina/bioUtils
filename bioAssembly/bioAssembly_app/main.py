#!/usr/bin/env python3
"""
Downstream analysis of assembled contigs.
Trim, map and assemble.
"""

# rich dependencies
from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich import box

import argparse, sys
#from .cli import final_colored_ascii_banner


parser = argparse.ArgumentParser(
    #description=final_colored_ascii_banner,
    add_help=False,  
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                    help='Show this help message and exit.')

# Required arguments
parser.add_argument('--inputdir', '-i', required=True)
parser.add_argument('--outdir', '-o', required=True)
parser.add_argument('--layout', required=True, choices=['single', 'paired'])


# Modules arguments 
parser.add_argument('--skip-qc', action='store_true')
parser.add_argument('--skip-trim', action='store_true')
parser.add_argument('--trim', choices=['fastp', 'trimmomatic'], required=True)
parser.add_argument('--skip-align', action='store_true')
parser.add_argument('--assembly', choices=['megahit', 'spades', 'rnaspades', 'trinity'], required=True)
parser.add_argument('--integrative-assembly', default='MT') 
# the initial letter of each assembly tool to be combined: 'MT', SRT', 'MSRT', etc.

# Temporary files - By default all the temp files are saved after processing
parser.add_argument('--remove-all', action='store_true')
parser.add_argument('--remove-trimmed', action='store_true')
parser.add_argument('--remove-aligned', action='store_true')
parser.add_argument('--remove-assembled', action='store_true')

# Config args
parser.add_argument('--threads', '-t', type=int, default=4)
parser.add_argument('--memory', '-m', type=int, default=8)


# rich console
console = Console()


def show_rich_help():
    """Display rich help message using Rich panels."""
    
    # You can customize your banner here if you have one
    # final_colored_ascii_banner = pyfiglet.figlet_format("YourToolName") 
    # print(final_colored_ascii_banner) 
    
    help_text = Text()
    help_text.append("\nA tool for sequence data processing and assembly.\n", style="italic")
    help_text.append("More info: https://github.com/gabrielvpina/bioutils\n", style="dim blue underline")
    help_text.append("\nQuality Control -> Alignment -> Assembly\n", style="bold green")
    console.print(help_text)

    # Required arguments section
    required_panel = Panel(
        "[bold white]--inputdir/-i[/bold white]\n"
        "  Input directory containing raw sequencing data.\n\n"
        "[bold white]--outdir/-o[/bold white]\n"
        "  Output directory to save all generated files.\n\n"
        "[bold white]--layout[/bold white]\n"
        "  If the libraries in input directory is single-end or paired-end.",
        title="[bold cyan]REQUIRED ARGUMENTS[/bold cyan]",
        border_style="cyan",
        width=70,
        box=box.ROUNDED
    )
    console.print(required_panel)
    console.print()

    # Module arguments section
    module_panel = Panel(
        "[bold white]--skip-qc[/bold white]\n"  
        "   Skip the quality control step.\n\n"
        "[bold white]--skip-trim[/bold white]\n"  
        "   Skip the read trimming step.\n\n"
        "[bold white]--trim[/bold white] [choices: fastp, trimmomatic]\n"  
        "   Specify the trimming tool to use [bold white](fastp or trimmomatic)[/bold white].\n\n"
        "[bold white]--skip-align[/bold white]\n" 
        "   Skip the alignment step.\n\n"
        "[bold white]--assembly[/bold white] \n"  
        "   Specify the assembly tool to use [bold white](megahit, spades, rnaspades, trinity)[/bold white].\n\n"
        "[bold white]--integrative-assembly[/bold white] [default: MT]\n"  
        "   Combine initial letters of assembly tools for integrative assembly (e.g., 'MT', 'SRT').",
        title="[bold green]MODULES ARGUMENTS[/bold green]",
        border_style="green",
        width=70,
        box=box.ROUNDED
    )
    console.print(module_panel)
    console.print()

    # Temporary files arguments section
    temp_files_panel = Panel(
        "[bold white]--remove-all[/bold white]\n"  
        "   Remove all intermediate temporary files after processing.\n\n"
        "[bold white]--remove-trimmed[/bold white]\n"  
        "   Remove trimmed reads after processing.\n\n"
        "[bold white]--remove-aligned[/bold white]\n"  
        "   Remove aligned reads after processing.\n\n"
        "[bold white]--remove-assembled[/bold white]\n"  
        "   Remove assembled contigs after processing.",
        title="[bold cyan]TEMPORARY FILES (OPTIONAL)[/bold cyan]",
        border_style="cyan",
        width=70,
        box=box.ROUNDED
    )
    console.print(temp_files_panel)
    console.print()

    # Config arguments section
    config_panel = Panel(
        "[bold white]--threads[/bold white] [type: int, default: 4]\n"  
        "   Number of CPU threads to use.\n\n"
        "[bold white]--memory[/bold white] [type: int, default: 8]\n"  
        "   Amount of memory in GB to allocate.",
        title="[bold green]CONFIGURATION ARGUMENTS[/bold green]",
        border_style="green",
        width=70,
        box=box.ROUNDED
    )
    console.print(config_panel)
    console.print()
    
    # Other options section (help, version)
    help_footer = Panel(
        "[bold white]--help / -h[/bold white]" 
        "   Show this help message and exit",
        title="[bold cyan]OTHER OPTIONS[/bold cyan]",
        border_style="cyan",
        width=70,
        box=box.ROUNDED
    )
    console.print(help_footer)

# Show message if no arguments are provided
if len(sys.argv) == 1:
    console.print(
        "[bold red]ERROR:[/bold red] No arguments provided. Use '-h' or '--help' to see available options.",
        style="red"
    )
    sys.exit(1)

# Check for help before creating parser
if "-h" in sys.argv or "--help" in sys.argv:
    show_rich_help()
    sys.exit(0)

# Parse normal
args = parser.parse_args()

if __name__ == "__main__":
    # Show message if no arguments are provided
    if len(sys.argv) == 1:
        console.print(
            "[bold red]ERROR:[/bold red] No arguments provided. Use '-h' or '--help' to see available options.",
            style="red"
        )
        sys.exit(1)
    # Check for help before parsing
    if "-h" in sys.argv or "--help" in sys.argv:
        show_rich_help(parser)
        sys.exit(0)

    # Create directory
    import os
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)



    # Parse arguments normally if help is not requested
    args = parser.parse_args()