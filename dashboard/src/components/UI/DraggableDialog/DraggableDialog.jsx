import React from 'react';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogContentText from '@material-ui/core/DialogContentText';
import DialogTitle from '@material-ui/core/DialogTitle';
import Paper from '@material-ui/core/Paper';
import Draggable from 'react-draggable';

import FlexContainer from "../../UI/FlexContainer/FlexContainer";


function PaperComponent(props) {
  return (
    <Draggable handle="#draggable-dialog-title" cancel={'[class*="MuiDialogContent-root"]'}>
      <Paper {...props} />
    </Draggable>
  );
}

const DraggableDialog = (props) => {
  const [index, setIndex] = React.useState(null);

  const handleItemClicked = (idx) => {
    setIndex(idx);
  }

  // console.log(props, "[DraggableDialog]")
  const content = <FlexContainer items={props.items} click={handleItemClicked} />
  return (
    <div key="dialog">
      <Dialog
        open={true}
        PaperComponent={PaperComponent}
        aria-labelledby="draggable-dialog-title"
      >
        <DialogTitle style={{ cursor: 'move' }} id="draggable-dialog-title">
          {props.title}
        </DialogTitle>
        <DialogContent>
          <DialogContentText>
            {content}
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button autoFocus onClick={ () => props.click(null) } color="primary">
            Cancel
          </Button>
          <Button onClick={ () => props.click(index) } color="primary">
            Accept
          </Button>
        </DialogActions>
      </Dialog>
    </div>
  );
}

export default DraggableDialog;