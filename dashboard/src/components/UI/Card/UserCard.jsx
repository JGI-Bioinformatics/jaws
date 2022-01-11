import React from 'react';

import Card from '@material-ui/core/Card';
import CardActions from '@material-ui/core/CardActions';
import CardContent from '@material-ui/core/CardContent';
import Button from '@material-ui/core/Button';
import Typography from '@material-ui/core/Typography';

import classes from "./UserCard.module.css"

const UserCard = props => {
  return (
    <Card className={classes.CardRoot}>
      <CardContent>
        <Typography variant="h5" component="h2">
          {props.name} : {props.email}
        </Typography>
        <Typography className={classes.CardPos} color="textSecondary">
          {props.id}
        </Typography>
        <Typography variant="body2" component="p">
          well meaning and kindly.
          <br />
          {'"a benevolent smile"'}
        </Typography>
      </CardContent>
      <CardActions>
        <Button size="small">Learn More</Button>
      </CardActions>
    </Card>
  )
}
export default UserCard;