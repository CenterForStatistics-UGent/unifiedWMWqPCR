# Internal function to show the group order
# Need a vector of 2, as contained in groupinfo slot

groups.message <- function(groups)
{
  group.left <- groups[1]
  group.right <- groups[2]
  prob.print <- paste("Fitted probabilities: P(",group.left," < ", group.right, ") + 0.5 P(", group.left," = ", group.right,")", sep ="")
  return(prob.print)
}