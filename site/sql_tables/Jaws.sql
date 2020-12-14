CREATE TABLE IF NOT EXISTS users (
  id CHAR(36) PRIMARY KEY,
  name CHAR(64),
  email CHAR(64) UNIQUE,
  is_admin BOOLEAN NOT NULL,
  auth_access_token CHAR(256),
  auth_refresh_token CHAR(256),
  auth_expires_at_seconds INTEGER,
  transfer_access_token CHAR(256),
  transfer_refresh_token CHAR(256),
  transfer_expires_at_seconds INTEGER,
  groups_access_token CHAR(256),
  groups_refresh_token CHAR(256),
  groups_expires_at_seconds INTEGER
);

CREATE TABLE IF NOT EXISTS workflows (
  id INTEGER PRIMARY KEY,
  name CHAR(32) NOT NULL,
  version CHAR(16) NOT NULL DEFAULT 'latest',
  user_id CHAR(36) NOT NULL,
  FOREIGN KEY(user_id) REFERENCES users(id),
  created DATETIME NOT NULL,
  updated DATETIME NOT NULL,
  is_released BOOLEAN NOT NULL DEFAULT FALSE,
  is_deprecated BOOLEAN NOT NULL DEFAULT FALSE,
  wdl TEXT NOT NULL,
  doc TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS sites (
  id INTEGER PRIMARY KEY,
  endpoint CHAR(36) NOT NULL,
  basepath CHAR(256) NOT NULL,
  staging CHAR(256) NOT NULL,
  max_ram_gb INTEGER NOT NULL,
  max_transfer_gb INTEGER NOT NULL,
  has_dev_shm BOOLEAN NOT NULL
);

CREATE TABLE IF NOT EXISTS runs (
  id INTEGER PRIMARY KEY,
  user_id CHAR(36),
  FOREIGN KEY (user_id) REFERENCES users(id),
  site_id CHAR(8),
  FOREIGN KEY (site_id) REFERENCES sites(id),
  submission_uuid CHAR(36),
  status CHAR(16),
  cromwell_id CHAR(36),
  submitted DATETIME NOT NULL,
  updated DATETIME NOT NULL,
  upload_task_id CHAR(36) NOT NULL,
  download_task_id CHAR(36) NOT NULL,
  dest_endpoint CHAR(36),
  dest_path CHAR(256)
);