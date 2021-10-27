import React from 'react';
import CircularProgress from '@material-ui/core/CircularProgress';
import { useQuery } from "react-query";
import ReactMarkdown from 'react-markdown';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogContentText from '@material-ui/core/DialogContentText';
import DialogTitle from '@material-ui/core/DialogTitle';

const fetchWDLReadme = async (url) => {
  const resp = await fetch(url,
  {  
    method: 'GET',
    headers: {
      'Authorization': `Bearer ${process.env.REACT_APP_MY_PWD}`,
    }
 });
  return await resp.json();
}

export default function ReadmeContent(props) {
  const {repo_id, repo_name} = props;
  const [open, setOpen] = React.useState(false);
  const [scroll, setScroll] = React.useState('paper');

  const handleClickOpen = (scrollType) => () => {
    setOpen(true);
    setScroll(scrollType);
  };

  const handleClose = () => {
    setOpen(false);
  };

  const descriptionElementRef = React.useRef(null);
  React.useEffect(() => {
    if (open) {
      const { current: descriptionElement } = descriptionElementRef;
      if (descriptionElement !== null) {
        descriptionElement.focus();
      }
    }
  }, [open]);

  let api_readme_url = `https://code.jgi.doe.gov/api/v4/projects/${repo_id}/repository/files/README%2Emd?ref=master`

  // USEFUL NOTE: when you want to pass in an argument to a function called 
  // by useQuery, wrap the function with a function.
  const { data, error, status } = useQuery(api_readme_url, () => fetchWDLReadme(api_readme_url) );

  if (error) return "An error has occurred: " + error.message;

  let readme_contents = null
  if (status === 'success') {
      // test if there is a readme
      if (!data.content) {
        readme_contents = 'empty'
      }else{
        readme_contents = new TextDecoder().decode(Buffer.from(data.content, 'base64'));
      }
  }

  return (
    <div>
      <Button onClick={handleClickOpen('paper')} style={{color: "blue"}}>readme</Button>
      <Dialog
        open={open}
        onClose={handleClose}
        fullWidth={true}
        maxWidth={'lg'}
        scroll={scroll}
        aria-labelledby="scroll-dialog-title"
        aria-describedby="scroll-dialog-description"
      >
        <DialogTitle id="scroll-dialog-title">{repo_name}</DialogTitle>
        <DialogContent dividers={scroll === 'paper'}>
          <DialogContentText
            id="scroll-dialog-description"
            ref={descriptionElementRef}
            tabIndex={-1}
            style={{color: "black"}}
          >
            {status === "loading" ? <CircularProgress /> : <ReactMarkdown>{readme_contents}</ReactMarkdown>}
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button onClick={handleClose} color="primary">
            Cancel
          </Button>
        </DialogActions>
      </Dialog>
    </div>
  );
}