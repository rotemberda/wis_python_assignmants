import tkinter as tk
from random import randint
from tkinter import font     # to set the font

# Constants
TEMP_LABEL_DURATION = 4000  # Milliseconds
WIN_MESSAGE_DURATION = 3000
MAX_SECRET_NUMBER = 20
CHEAT_DURATION = 2000  # 2 seconds for cheat message

# scoreboard class
class GameStats:
    def __init__(self):
        self.rounds_won = 0
        self.total_guesses = 0
    
    @property
    def average_guesses(self):
        if self.rounds_won == 0:
            return 0
        return round(self.total_guesses / self.rounds_won, 1)
    
    def update(self, guesses):
        """update the game stats"""
        self.rounds_won += 1
        self.total_guesses += guesses

    def penalize(self, penalty):
        """Add penalty guesses to total without incrementing rounds won"""
        self.total_guesses += penalty



'''
Create the main application window
'''
# Initialize the app
app = tk.Tk()
app.title("Guessing Game")

# Initialize game statistics
stats = GameStats()

# Maximize game window om startup
app.state('zoomed')   # For Windows
# OR
# app.attributes('-zoomed', True)  # For Linux

app.resizable(True, True) # enable the user to resize the window



'''
Aesthetics
'''
# Set a larger global font
default_font = font.nametofont("TkDefaultFont")
default_font.configure(size=16)  # Set font size to 16
app.option_add("*Font", default_font)  # Apply the font globally

# Create control buttons frame (top-left)
control_frame = tk.Frame(app)
control_frame.pack(side="top", anchor="nw", padx=10, pady=10)

# Create scoreboard frame (bottom-left)
scoreboard_frame = tk.Frame(app)
scoreboard_frame.pack(side="bottom", anchor="sw", padx=10, pady=10)

rounds_label = tk.Label(scoreboard_frame, text="Rounds Won: 0", font=("TkDefaultFont", 12))
rounds_label.pack(anchor="w")

average_label = tk.Label(scoreboard_frame, text="Average Guesses: 0", font=("TkDefaultFont", 12))
average_label.pack(anchor="w")



'''
Useful functions
'''
def update_scoreboard():
    """Update the scoreboard labels with current stats"""

    rounds_label.config(text=f"Rounds Won: {stats.rounds_won}")
    average_label.config(text=f"Average Guesses: {stats.average_guesses}")


def quit_game():
    """Clean screen and display quit message"""

    for widget in app.winfo_children():
        widget.destroy()
    
    quit_label = tk.Label(app, text="That was fun!!\n\nQuiting game...")
    quit_label.pack(expand=True)
    app.after(2000, app.quit)  # Quit after 2 seconds


def temporary_label(parent, text, duration=TEMP_LABEL_DURATION):
    """Display a temporary label with the specified text for a limited duration"""

    temp_label = tk.Label(parent, text=text)
    temp_label.pack()
    parent.after(duration, temp_label.pack_forget)



'''
Game logic
'''
def game_round():
    """Start a new game round."""   

    # Hide the welcoming prompt and "Start a Round" button if they exist
    for widget in app.winfo_children():
        if widget not in [control_frame, scoreboard_frame]:
            widget.pack_forget()

    # Create round-specific data dictionary
    round_data = {
        'secret_number': randint(1, MAX_SECRET_NUMBER),  # Generate a random number
        'guess_num': 0  # Initialize guess counter
    }

    # Create frame and widgets
    button_frame = tk.Frame(app)
    button_frame.pack(pady=10)

    prompt_guess = tk.Label(app, text=f'Guess a number from 1 to {MAX_SECRET_NUMBER}:')
    prompt_guess.pack()

    # Store widgets in round_data
    round_data.update({
        'button_frame': button_frame,
        'prompt_guess': prompt_guess
    })

    """Create control buttons for the round"""
    quit_round_btn = tk.Button(control_frame, text="Quit Round",
                              command=lambda: quit_round(round_data),
                              font=("TkDefaultFont", 12))
    quit_round_btn.pack(side="top", anchor="w", pady=2)
    
    cheat_btn = tk.Button(control_frame, text="Cheat",
                         command=lambda: show_secret(round_data),
                         font=("TkDefaultFont", 12))
    cheat_btn.pack(side="top", anchor="w", pady=2)

    """Buttons for guess"""
    for i in range(1, MAX_SECRET_NUMBER + 1):
        if i == round_data['secret_number']:
            btn = tk.Button(button_frame, text=str(i),
                          command=lambda: right_guess(round_data),
                          activebackground="green", activeforeground="white")
        else:
            btn = tk.Button(button_frame, text=str(i),
                          command=wrong_guess(i, round_data),
                          activebackground="red", activeforeground="white")
        btn.pack(side="left", padx=2)



'''
Helper functions
'''
def right_guess(round_data):
    """Handle correct guess"""

    round_data['guess_num'] += 1
    round_data['button_frame'].destroy()  # Remove all buttons
    stats.update(round_data['guess_num'])  # Update game statistics
    update_scoreboard()  # Update the scoreboard display
    display_win_message(round_data)


def wrong_guess(guess_value, round_data):
    """Create handler for incorrect guesses"""

    def handler():
        round_data['guess_num'] += 1
        if guess_value < round_data['secret_number']:
            temporary_label(app, "Your guess is too small! 😔 Try again.")
        elif guess_value > round_data['secret_number']:
            temporary_label(app, "Your guess is too big! 😔 Try again.")
    return handler


def start_new_round():
    """Prompt for a new game round."""

    for widget in control_frame.winfo_children():
        if widget != quit_game_btn:  # Keep only the quit game button
            widget.destroy()

    round_prompt = tk.Label(app, text="Let's have another round?")  # Prompt for another round
    round_prompt.pack()


    def restart_game():
        """Clean the screen and restarting a new game."""

        round_prompt.pack_forget()
        new_round_button.pack_forget()
        game_round()

    new_round_button = tk.Button(
        app, text="Start a New Round", command=restart_game,
        activebackground="blue", activeforeground="white"
    )
    new_round_button.pack()


def display_win_message(round_data):
    """Display the winning message"""

    round_data['prompt_guess'].pack_forget()  # Hide current round widgets
    winning_message = tk.Label(app, text=f"You Won!😊\n\nIt took you {round_data['guess_num']} guesses.")
    winning_message.pack()
    app.after(WIN_MESSAGE_DURATION, lambda: [winning_message.pack_forget(), start_new_round()])



'''
Handeling the controle buttons
'''
def quit_round(round_data):
    """Handle quitting the current round"""

    round_data['button_frame'].destroy()
    round_data['prompt_guess'].pack_forget()
    stats.penalize(2)  # Add 2 penalty guesses
    update_scoreboard()
    
    # Display message and start new round
    quit_round_label = tk.Label(app, text="Round ended early!\nPenalty: +2 guesses")
    quit_round_label.pack()
    app.after(2000, lambda: [quit_round_label.pack_forget(), start_new_round()])


def show_secret(round_data):
    """Show the secret number temporarily"""

    temporary_label(app, f"The secret number is {round_data['secret_number']}, shhh...", CHEAT_DURATION)



'''
Initializing game
'''
# Create initial control buttons
quit_game_btn = tk.Button(control_frame, text="Quit Game", command=quit_game,
                         font=("TkDefaultFont", 12))
quit_game_btn.pack(side="top", anchor="w", pady=2)

# Initial widgets
welcoming_prompt = tk.Label(app, text="Welcome to the guessing game!!")
welcoming_prompt.pack()

start_round = tk.Button(
    app, text="Start a Round", width=25, command=game_round,
    activebackground="blue", activeforeground="white"
)
start_round.pack()

app.mainloop()
