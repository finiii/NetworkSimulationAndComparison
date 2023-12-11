#this is for evaluating the results from coda4microbiome

count <- 0

# Iterate through the list and count values greater than 0.05
for (i in c(1:100)) {
  number <- temp[[i]]$p_value_coda
  if (number < 0.05) {
    count <- count + 1
  }
}

# Display the count
print(count)