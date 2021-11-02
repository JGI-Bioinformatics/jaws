import React from 'react';

import { makeStyles } from '@material-ui/core/styles';
import Pagination from '@material-ui/lab/Pagination';

const useStyles = makeStyles((theme) => ({
  root: {
    display: "flex",
    '& > *': {
      marginTop: theme.spacing(2),
    },
  },
  pager: {
    marginLeft: "auto"
  }
}));

const MyPagination = props => {
  const classes = useStyles();
  
  const { initpage, onchange, count} = props;

  const [page, setPage] = React.useState(+initpage);

  const handleChangePage = (event, newPage) => {
    setPage(newPage);
    onchange(newPage);
  };

  return (
    <div className={classes.root}>
        <Pagination count={+count} 
                    page={+page}
                    onChange={handleChangePage}
                    variant="outlined" 
                    color="secondary" 
                    className={classes.pager}/>
      </div>
  )
}

export default MyPagination;