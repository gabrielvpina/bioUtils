

# this is pure LLM code 
def get_gradient_color(progress, colors):
    """
    Interpolates a color along a gradient defined by a list of RGB colors.
    """
    if progress <= 0:
        return colors[0]
    if progress >= 1:
        return colors[-1]

    num_segments = len(colors) - 1
    segment_index = int(progress * num_segments)
    
    if segment_index == num_segments:
        segment_index -= 1 

    start_color = colors[segment_index]
    end_color = colors[segment_index + 1]

    # Calculate progress within the current segment
    segment_progress = (progress * num_segments) - segment_index

    # Interpolate each RGB component
    r = int(start_color[0] + segment_progress * (end_color[0] - start_color[0]))
    g = int(start_color[1] + segment_progress * (end_color[1] - start_color[1]))
    b = int(start_color[2] + segment_progress * (end_color[2] - start_color[2]))
    
    return (r, g, b)


# The ASCII art banner provided by the user
ascii_banner = f""" 
▗▄▄▖ ▗▄▄▄▖ ▗▄▖  ▗▄▖  ▗▄▄▖ ▗▄▄▖▗▄▄▄▖▗▖  ▗▖▗▄▄▖ ▗▖ ▗▖  ▗▖
▐▌ ▐▌  █  ▐▌ ▐▌▐▌ ▐▌▐▌   ▐▌   ▐▌   ▐▛▚▞▜▌▐▌ ▐▌▐▌  ▝▚▞▘ 
▐▛▀▚▖  █  ▐▌ ▐▌▐▛▀▜▌ ▝▀▚▖ ▝▀▚▖▐▛▀▀▘▐▌  ▐▌▐▛▀▚▖▐▌   ▐▌  
▐▙▄▞▘▗▄█▄▖▝▚▄▞▘▐▌ ▐▌▗▄▄▞▘▗▄▄▞▘▐▙▄▄▖▐▌  ▐▌▐▙▄▞▘▐▙▄▄▖▐▌  
"""



# Split the banner into individual lines
lines = ascii_banner.strip().split('\n')

# Determine the maximum width of the banner to properly apply the gradient
max_width = max(len(line) for line in lines)

# Define the gradient colors (RGB tuples)
# These colors are chosen to approximate the blue -> purple -> pink/red gradient from the image.
gradient_colors = [
    (144, 238, 144),  # Verde Claro (Light Green)
    (60, 179, 113),   # Verde Médio (Medium Green)
    (34, 139, 34)     # Verde Escuro (Forest Green)
]


# List to store the fully colored lines
colored_banner_lines = []

# Iterate through each line and each character to apply the gradient
for line in lines:
    colored_line_chars = []
    for i, char in enumerate(line):
        # Do not color spaces; they should remain transparent to show the terminal background
        if char == ' ': 
            colored_line_chars.append(' ')
            continue
        
        # Calculate the progress for the current character across the entire banner width
        # This determines where in the gradient the character's color should be
        progress = i / (max_width - 1) if max_width > 1 else 0
        
        # get the interpolated RGB color for the current progress
        r, g, b = get_gradient_color(progress, gradient_colors)
        
        # Apply ANSI 24-bit true color escape code to the character
        # \x1b[38;2;R;G;Bm sets the foreground color
        colored_line_chars.append(f"\x1b[38;2;{r};{g};{b}m{char}")
    
    # reset color at the end of each line to avoid bleeding into subsequent terminal output
    # \x1b[0m resets all SGR (Select Graphic Rendition) attributes
    colored_banner_lines.append("".join(colored_line_chars) + "\x1b[0m") 

final_colored_ascii_banner = "\n".join(colored_banner_lines)
#final_colored_ascii_banner = print("\n".join(colored_banner_lines))
# print(final_colored_ascii_banner)