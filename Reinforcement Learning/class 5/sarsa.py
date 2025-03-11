import gymnasium as gym
import numpy as np
from matplotlib import pyplot as plt
import time

env = gym.make('CliffWalking-v0')

q_table = np.zeros((48, 4))

def policy(state, explore=0.0):
    action = int(np.argmax(q_table[state]))
    if np.random.random() <= explore:
        action = int(np.random.randint(low=0, high=4, size=1))
    return action


# Parameters
EPSILON = 0.1
ALPHA = 0.1
GAMMA = 0.9

NUM_EPISODE = 500

episodes_rewards = []

for episode in range(NUM_EPISODE):

    done = False

    total_reward = 0
    episode_length = 0

    state = env.reset()[0]
    action = policy(state, EPSILON)

    while not done:
        step = env.step(action)
        next_state, reward, done, _ = step[0], step[1], step[2], step[3]
        next_action = policy(next_state, EPSILON)

        q_table[state][action] += ALPHA * (reward + GAMMA * q_table[next_state][next_action]-q_table[state][action])

        state = next_state
        action = next_action

        total_reward += reward
        episode_length += 1
    
    episodes_rewards.append(total_reward)

    print("Episode:", episode, "Episode Length: ",
          episode_length, "Total Reward: ", total_reward)

env.close()

print("Training finished !")

# Show the result of our training
EPSILON = 0.0
env = gym.make('CliffWalking-v0', render_mode="human")
done = False

state = env.reset()[0]
action = policy(state, EPSILON)

while not done:
    time.sleep(1)

    step = env.step(action)
    next_state, reward, done, _ = step[0], step[1], step[2], step[3]
    next_action = policy(next_state, EPSILON)
    state = next_state
    action = next_action

    total_reward += reward
    episode_length += 1

plt.plot([i for i in range(len(episodes_rewards))], episodes_rewards)
plt.show()