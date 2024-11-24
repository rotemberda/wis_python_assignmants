from random import randint
import sys


def main():
    # print a message to the user about the start of the game
    print("Welcome to the guessing game!!")
    
    
    while True:

        # start a round
        game_round()

        #round ended

        # ask if user wants another round
        another_round = input("That was fun:) \nDo you want to play another round? (Yes/No) ").lower()

        if another_round == 'yes' or another_round == 'y':
            continue
        
        # user doesn't want anopther round. quiting the game with a nice message
        else:
            print('\n' * 10, "See you next time :)", 'Quiting game...', sep = "\n")
            break




def game_round():               # function that starts a game round

    # generate a random number
    global secret_number
    secret_number = randint(1, 20)

    #setting the guesses number to zero
    guess_num = 0
   
    while True:
        # ask for the user's guess
        guess = prompt_guess()
        
        # check if user wanted to quite the round
        if guess == 'n':
            print("Ending round.")
            break
    
        # compare to the secret number and inform the user
        if guess == secret_number:
            print("You Won!", f"You needed {guess_num} guesses.", sep = "\n" * 2)
            break

        # the user's guess wasn't right
        # checking to see if it was too small or too big, adding 1 to the guesses number
        else:   
            print(check_value(guess))
            guess_num += 1

    


def prompt_guess():            # ask for user's guess and check if valid
    
    while True:

        guess = input("What's your guess? ")
        
    # check for special inputs

        if guess.lower() == 'x':
            print('\n' * 10, "See you next time :)", 'Quiting game...', sep = "\n")
            sys.exit(0)

        if guess == 'n':
            return guess.lower()
        
        elif guess.lower() == 's':
            print("I usually do not condone cheating, but it seems this is a tough one...", f"the secrete number is {secret_number}", sep = "\n" * 2, end = '\n' * 2)
            return prompt_guess()
        
    # validating the user's input

        try:
            return int(guess)
        
        except:
            print("Invalid input! Input type isn't 'int'")
    

    


def check_value(guess):     # check whether the user's guess is too small or too big    
    
    if guess < secret_number:
        return "Your guess is too small! Try again."
    
    else:
        return "Your guess is too big! Try again."
    



if __name__ == '__main__':
    main()  